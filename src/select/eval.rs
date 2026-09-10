//! Compilation of the AST against a column map and per-record evaluation.

use std::collections::HashMap;

use crate::error::{KiraError, Result};
use crate::pairs::columns::ColumnMap;
use crate::pairs::record::PairRecordRef;
use crate::select::parser::{Ast, BinOp, parse};
use crate::util::int::{parse_f64, parse_i64};

/// Column value type used when reading a column.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ColumnType {
    /// Parsed as a signed integer.
    Int,
    /// Parsed as a float.
    Float,
    /// Raw bytes.
    Str,
}

/// Runtime value.
#[derive(Debug, Clone, PartialEq)]
pub enum Value<'a> {
    /// Integer.
    Int(i64),
    /// Float.
    Float(f64),
    /// Borrowed string bytes.
    Str(&'a [u8]),
    /// Owned string bytes (function results).
    OwnedStr(Vec<u8>),
    /// Boolean.
    Bool(bool),
    /// List (from literals).
    List(Vec<Value<'a>>),
}

impl Value<'_> {
    fn truthy(&self) -> bool {
        match self {
            Value::Int(i) => *i != 0,
            Value::Float(f) => *f != 0.0,
            Value::Str(s) => !s.is_empty(),
            Value::OwnedStr(s) => !s.is_empty(),
            Value::Bool(b) => *b,
            Value::List(l) => !l.is_empty(),
        }
    }

    fn as_bytes(&self) -> Option<&[u8]> {
        match self {
            Value::Str(s) => Some(s),
            Value::OwnedStr(s) => Some(s),
            _ => None,
        }
    }

    fn type_name(&self) -> &'static str {
        match self {
            Value::Int(_) => "int",
            Value::Float(_) => "float",
            Value::Str(_) | Value::OwnedStr(_) => "str",
            Value::Bool(_) => "bool",
            Value::List(_) => "list",
        }
    }
}

/// Compiled expression node.
#[derive(Debug, Clone)]
enum Expr {
    Int(i64),
    Float(f64),
    Str(Vec<u8>),
    Bool(bool),
    Column(usize, ColumnType, String),
    List(Vec<Expr>),
    Neg(Box<Expr>),
    Not(Box<Expr>),
    And(Box<Expr>, Box<Expr>),
    Or(Box<Expr>, Box<Expr>),
    Bin(BinOp, Box<Expr>, Box<Expr>),
    Abs(Box<Expr>),
    ToInt(Box<Expr>),
    ToFloat(Box<Expr>),
    ToStr(Box<Expr>),
    Len(Box<Expr>),
    Min(Vec<Expr>),
    Max(Vec<Expr>),
    CsvMatch(Box<Expr>, Vec<Vec<u8>>),
    WildcardMatch(Box<Expr>, regex::bytes::Regex),
    RegexMatch(Box<Expr>, regex::bytes::Regex),
}

/// A compiled filter.
#[derive(Debug, Clone)]
pub struct Filter {
    expr: Expr,
    max_col: usize,
}

fn default_types() -> HashMap<String, ColumnType> {
    let mut m = HashMap::new();
    for c in ["pos1", "pos2", "mapq1", "mapq2"] {
        m.insert(c.to_string(), ColumnType::Int);
    }
    m
}

/// Parse a `--type-cast` type name.
pub fn parse_column_type(s: &str) -> Result<ColumnType> {
    match s {
        "int" => Ok(ColumnType::Int),
        "float" => Ok(ColumnType::Float),
        "str" => Ok(ColumnType::Str),
        other => Err(KiraError::arg(format!(
            "unknown type {other:?} (expected int, float or str)"
        ))),
    }
}

fn glob_to_regex(glob: &str) -> String {
    let mut re = String::from("^");
    for c in glob.chars() {
        match c {
            '*' => re.push_str(".*"),
            '?' => re.push('.'),
            '[' => re.push('['),
            ']' => re.push(']'),
            other => re.push_str(&regex::escape(&other.to_string())),
        }
    }
    re.push('$');
    re
}

impl Filter {
    /// Compile an expression for the given columns. `type_cast` overrides
    /// the default column types (`pos1`, `pos2`, `mapq1`, `mapq2` are int).
    pub fn compile(
        src: &str,
        cols: &ColumnMap,
        type_cast: &[(String, ColumnType)],
    ) -> Result<Self> {
        let ast = parse(src)?;
        let mut types = default_types();
        for (c, t) in type_cast {
            types.insert(c.clone(), *t);
        }
        let mut max_col = 0usize;
        let expr = compile(&ast, cols, &types, &mut max_col)?;
        Ok(Self { expr, max_col })
    }

    /// Highest column index referenced.
    pub fn max_column(&self) -> usize {
        self.max_col
    }

    /// Evaluate against a record.
    #[inline]
    pub fn eval(&self, rec: &PairRecordRef<'_>) -> Result<bool> {
        Ok(eval(&self.expr, rec)?.truthy())
    }
}

fn compile(
    ast: &Ast,
    cols: &ColumnMap,
    types: &HashMap<String, ColumnType>,
    max_col: &mut usize,
) -> Result<Expr> {
    let mut sub = |a: &Ast| compile(a, cols, types, max_col);
    Ok(match ast {
        Ast::Int(i) => Expr::Int(*i),
        Ast::Float(f) => Expr::Float(*f),
        Ast::Str(s) => Expr::Str(s.as_bytes().to_vec()),
        Ast::Bool(b) => Expr::Bool(*b),
        Ast::Name(n) => match cols.names().iter().position(|c| c == n) {
            Some(i) => {
                *max_col = (*max_col).max(i);
                Expr::Column(
                    i,
                    types.get(n).copied().unwrap_or(ColumnType::Str),
                    n.clone(),
                )
            }
            None => {
                return Err(KiraError::Expression(format!(
                    "unknown column {n:?} (known columns: {})",
                    cols.names().join(", ")
                )));
            }
        },
        Ast::ColIndex(i) => {
            *max_col = (*max_col).max(*i);
            let name = cols.names().get(*i).cloned().unwrap_or_default();
            Expr::Column(
                *i,
                types.get(&name).copied().unwrap_or(ColumnType::Str),
                format!("COLS[{i}]"),
            )
        }
        Ast::List(items) => Expr::List(items.iter().map(&mut sub).collect::<Result<_>>()?),
        Ast::Neg(e) => Expr::Neg(Box::new(sub(e)?)),
        Ast::Not(e) => Expr::Not(Box::new(sub(e)?)),
        Ast::And(a, b) => Expr::And(Box::new(sub(a)?), Box::new(sub(b)?)),
        Ast::Or(a, b) => Expr::Or(Box::new(sub(a)?), Box::new(sub(b)?)),
        Ast::Bin(op, a, b) => Expr::Bin(*op, Box::new(sub(a)?), Box::new(sub(b)?)),
        Ast::Call(name, args) => {
            let arity = |n: usize| -> Result<()> {
                if args.len() != n {
                    Err(KiraError::Expression(format!(
                        "{name}() takes {n} argument(s), got {}",
                        args.len()
                    )))
                } else {
                    Ok(())
                }
            };
            let literal = |a: &Ast| -> Result<String> {
                match a {
                    Ast::Str(s) => Ok(s.clone()),
                    _ => Err(KiraError::Expression(format!(
                        "{name}() requires a string literal as its second argument"
                    ))),
                }
            };
            match name.as_str() {
                "abs" => {
                    arity(1)?;
                    Expr::Abs(Box::new(sub(&args[0])?))
                }
                "int" => {
                    arity(1)?;
                    Expr::ToInt(Box::new(sub(&args[0])?))
                }
                "float" => {
                    arity(1)?;
                    Expr::ToFloat(Box::new(sub(&args[0])?))
                }
                "str" => {
                    arity(1)?;
                    Expr::ToStr(Box::new(sub(&args[0])?))
                }
                "len" => {
                    arity(1)?;
                    Expr::Len(Box::new(sub(&args[0])?))
                }
                "min" | "max" => {
                    if args.is_empty() {
                        return Err(KiraError::Expression(format!(
                            "{name}() needs at least one argument"
                        )));
                    }
                    let items = args.iter().map(&mut sub).collect::<Result<Vec<_>>>()?;
                    if name == "min" {
                        Expr::Min(items)
                    } else {
                        Expr::Max(items)
                    }
                }
                "csv_match" => {
                    arity(2)?;
                    let csv = literal(&args[1])?;
                    Expr::CsvMatch(
                        Box::new(sub(&args[0])?),
                        csv.split(',').map(|s| s.as_bytes().to_vec()).collect(),
                    )
                }
                "wildcard_match" => {
                    arity(2)?;
                    let pat = glob_to_regex(&literal(&args[1])?);
                    let re = regex::bytes::Regex::new(&pat)
                        .map_err(|e| KiraError::Expression(format!("bad wildcard: {e}")))?;
                    Expr::WildcardMatch(Box::new(sub(&args[0])?), re)
                }
                "regex_match" => {
                    arity(2)?;
                    let pat = format!("^(?:{})$", literal(&args[1])?);
                    let re = regex::bytes::Regex::new(&pat)
                        .map_err(|e| KiraError::Expression(format!("bad regex: {e}")))?;
                    Expr::RegexMatch(Box::new(sub(&args[0])?), re)
                }
                other => {
                    return Err(KiraError::Expression(format!(
                        "unknown function {other}() (supported: abs, int, float, str, len, min, max, csv_match, wildcard_match, regex_match)"
                    )));
                }
            }
        }
    })
}

fn type_error(what: &str, a: &Value<'_>, b: &Value<'_>) -> KiraError {
    KiraError::Expression(format!(
        "cannot apply {what} to {} and {} (hint: use --type-cast COLUMN int to compare numerically)",
        a.type_name(),
        b.type_name()
    ))
}

fn read_column<'a>(
    rec: &PairRecordRef<'a>,
    i: usize,
    t: ColumnType,
    name: &str,
) -> Result<Value<'a>> {
    let f = rec.field(i).ok_or_else(|| {
        KiraError::Expression(format!(
            "record has only {} columns, cannot read {name}",
            rec.n_fields()
        ))
    })?;
    Ok(match t {
        ColumnType::Str => Value::Str(f),
        ColumnType::Int => Value::Int(parse_i64(f).ok_or_else(|| {
            KiraError::Expression(format!(
                "column {name} value {:?} is not an integer",
                String::from_utf8_lossy(f)
            ))
        })?),
        ColumnType::Float => Value::Float(parse_f64(f).ok_or_else(|| {
            KiraError::Expression(format!(
                "column {name} value {:?} is not a number",
                String::from_utf8_lossy(f)
            ))
        })?),
    })
}

fn as_f64(v: &Value<'_>) -> Option<f64> {
    match v {
        Value::Int(i) => Some(*i as f64),
        Value::Float(f) => Some(*f),
        Value::Bool(b) => Some(if *b { 1.0 } else { 0.0 }),
        _ => None,
    }
}

fn as_i64(v: &Value<'_>) -> Option<i64> {
    match v {
        Value::Int(i) => Some(*i),
        Value::Bool(b) => Some(i64::from(*b)),
        _ => None,
    }
}

fn compare(op: BinOp, a: &Value<'_>, b: &Value<'_>) -> Result<bool> {
    use std::cmp::Ordering;
    let ord: Option<Ordering> = match (a, b) {
        (Value::Int(x), Value::Int(y)) => Some(x.cmp(y)),
        (Value::Bool(_), _) | (_, Value::Bool(_)) if matches!(op, BinOp::Eq | BinOp::Ne) => {
            let (x, y) = (as_f64(a), as_f64(b));
            match (x, y) {
                (Some(x), Some(y)) => x.partial_cmp(&y),
                _ => {
                    return Ok(!matches!(op, BinOp::Eq));
                }
            }
        }
        _ => {
            if let (Some(x), Some(y)) = (a.as_bytes(), b.as_bytes()) {
                Some(x.cmp(y))
            } else if let (Some(x), Some(y)) = (as_f64(a), as_f64(b)) {
                x.partial_cmp(&y)
            } else if matches!(op, BinOp::Eq | BinOp::Ne) {
                // Mixed str/number equality is simply false in Python.
                return Ok(op == BinOp::Ne);
            } else {
                return Err(type_error("comparison", a, b));
            }
        }
    };
    let Some(ord) = ord else {
        // NaN comparisons are false (except !=).
        return Ok(op == BinOp::Ne);
    };
    Ok(match op {
        BinOp::Eq => ord == Ordering::Equal,
        BinOp::Ne => ord != Ordering::Equal,
        BinOp::Lt => ord == Ordering::Less,
        BinOp::Le => ord != Ordering::Greater,
        BinOp::Gt => ord == Ordering::Greater,
        BinOp::Ge => ord != Ordering::Less,
        _ => false,
    })
}

fn arith<'a>(op: BinOp, a: Value<'a>, b: Value<'a>) -> Result<Value<'a>> {
    if let (Some(x), Some(y)) = (as_i64(&a), as_i64(&b)) {
        return Ok(match op {
            BinOp::Add => Value::Int(
                x.checked_add(y)
                    .ok_or_else(|| KiraError::Expression("integer overflow".into()))?,
            ),
            BinOp::Sub => Value::Int(
                x.checked_sub(y)
                    .ok_or_else(|| KiraError::Expression("integer overflow".into()))?,
            ),
            BinOp::Mul => Value::Int(
                x.checked_mul(y)
                    .ok_or_else(|| KiraError::Expression("integer overflow".into()))?,
            ),
            BinOp::Div => {
                if y == 0 {
                    return Err(KiraError::Expression("division by zero".into()));
                }
                Value::Float(x as f64 / y as f64)
            }
            BinOp::FloorDiv => {
                if y == 0 {
                    return Err(KiraError::Expression("division by zero".into()));
                }
                let q = x / y;
                let r = x % y;
                Value::Int(if r != 0 && ((r < 0) != (y < 0)) {
                    q - 1
                } else {
                    q
                })
            }
            BinOp::Mod => {
                if y == 0 {
                    return Err(KiraError::Expression("division by zero".into()));
                }
                let r = x % y;
                Value::Int(if r != 0 && ((r < 0) != (y < 0)) {
                    r + y
                } else {
                    r
                })
            }
            BinOp::Pow => {
                if y >= 0 {
                    Value::Int(
                        x.checked_pow(y as u32)
                            .ok_or_else(|| KiraError::Expression("integer overflow".into()))?,
                    )
                } else {
                    Value::Float((x as f64).powf(y as f64))
                }
            }
            _ => unreachable!(),
        });
    }
    if let (Value::Str(_) | Value::OwnedStr(_), Value::Str(_) | Value::OwnedStr(_)) = (&a, &b)
        && op == BinOp::Add
    {
        let mut v = a.as_bytes().unwrap_or(b"").to_vec();
        v.extend_from_slice(b.as_bytes().unwrap_or(b""));
        return Ok(Value::OwnedStr(v));
    }
    let (Some(x), Some(y)) = (as_f64(&a), as_f64(&b)) else {
        return Err(type_error("arithmetic", &a, &b));
    };
    Ok(Value::Float(match op {
        BinOp::Add => x + y,
        BinOp::Sub => x - y,
        BinOp::Mul => x * y,
        BinOp::Div => {
            if y == 0.0 {
                return Err(KiraError::Expression("division by zero".into()));
            }
            x / y
        }
        BinOp::FloorDiv => {
            if y == 0.0 {
                return Err(KiraError::Expression("division by zero".into()));
            }
            (x / y).floor()
        }
        BinOp::Mod => {
            if y == 0.0 {
                return Err(KiraError::Expression("division by zero".into()));
            }
            x - y * (x / y).floor()
        }
        BinOp::Pow => x.powf(y),
        _ => unreachable!(),
    }))
}

fn contains(item: &Value<'_>, list: &Value<'_>) -> Result<bool> {
    match list {
        Value::List(items) => {
            for it in items {
                if compare(BinOp::Eq, item, it)? {
                    return Ok(true);
                }
            }
            Ok(false)
        }
        Value::Str(_) | Value::OwnedStr(_) => {
            let hay = list.as_bytes().unwrap_or(b"");
            let needle = item.as_bytes().ok_or_else(|| {
                KiraError::Expression("'in <str>' requires a string on the left".into())
            })?;
            Ok(memchr::memmem::find(hay, needle).is_some())
        }
        _ => Err(KiraError::Expression(
            "'in' requires a list or string on the right".into(),
        )),
    }
}

fn eval<'a>(e: &Expr, rec: &PairRecordRef<'a>) -> Result<Value<'a>> {
    Ok(match e {
        Expr::Int(i) => Value::Int(*i),
        Expr::Float(f) => Value::Float(*f),
        Expr::Str(s) => Value::OwnedStr(s.clone()),
        Expr::Bool(b) => Value::Bool(*b),
        Expr::Column(i, t, name) => read_column(rec, *i, *t, name)?,
        Expr::List(items) => {
            Value::List(items.iter().map(|i| eval(i, rec)).collect::<Result<_>>()?)
        }
        Expr::Neg(x) => match eval(x, rec)? {
            Value::Int(i) => Value::Int(-i),
            Value::Float(f) => Value::Float(-f),
            Value::Bool(b) => Value::Int(-i64::from(b)),
            other => {
                return Err(KiraError::Expression(format!(
                    "cannot negate {}",
                    other.type_name()
                )));
            }
        },
        Expr::Not(x) => Value::Bool(!eval(x, rec)?.truthy()),
        Expr::And(a, b) => {
            let l = eval(a, rec)?;
            if !l.truthy() { l } else { eval(b, rec)? }
        }
        Expr::Or(a, b) => {
            let l = eval(a, rec)?;
            if l.truthy() { l } else { eval(b, rec)? }
        }
        Expr::Bin(op, a, b) => {
            let l = eval(a, rec)?;
            let r = eval(b, rec)?;
            match op {
                BinOp::Eq | BinOp::Ne | BinOp::Lt | BinOp::Le | BinOp::Gt | BinOp::Ge => {
                    Value::Bool(compare(*op, &l, &r)?)
                }
                BinOp::In => Value::Bool(contains(&l, &r)?),
                BinOp::NotIn => Value::Bool(!contains(&l, &r)?),
                _ => arith(*op, l, r)?,
            }
        }
        Expr::Abs(x) => match eval(x, rec)? {
            Value::Int(i) => Value::Int(i.abs()),
            Value::Float(f) => Value::Float(f.abs()),
            Value::Bool(b) => Value::Int(i64::from(b)),
            other => {
                return Err(KiraError::Expression(format!(
                    "abs() of {}",
                    other.type_name()
                )));
            }
        },
        Expr::ToInt(x) => match eval(x, rec)? {
            Value::Int(i) => Value::Int(i),
            Value::Float(f) => Value::Int(f.trunc() as i64),
            Value::Bool(b) => Value::Int(i64::from(b)),
            other => {
                let s = other.as_bytes().unwrap_or(b"");
                Value::Int(parse_i64(s).ok_or_else(|| {
                    KiraError::Expression(format!(
                        "int() cannot parse {:?}",
                        String::from_utf8_lossy(s)
                    ))
                })?)
            }
        },
        Expr::ToFloat(x) => match eval(x, rec)? {
            Value::Int(i) => Value::Float(i as f64),
            Value::Float(f) => Value::Float(f),
            Value::Bool(b) => Value::Float(if b { 1.0 } else { 0.0 }),
            other => {
                let s = other.as_bytes().unwrap_or(b"");
                Value::Float(parse_f64(s).ok_or_else(|| {
                    KiraError::Expression(format!(
                        "float() cannot parse {:?}",
                        String::from_utf8_lossy(s)
                    ))
                })?)
            }
        },
        Expr::ToStr(x) => match eval(x, rec)? {
            Value::Int(i) => Value::OwnedStr(i.to_string().into_bytes()),
            Value::Float(f) => {
                Value::OwnedStr(crate::util::pyfloat::format_py_float(f).into_bytes())
            }
            Value::Bool(b) => Value::OwnedStr(if b {
                b"True".to_vec()
            } else {
                b"False".to_vec()
            }),
            Value::Str(s) => Value::Str(s),
            Value::OwnedStr(s) => Value::OwnedStr(s),
            Value::List(_) => {
                return Err(KiraError::Expression(
                    "str() of a list is not supported".into(),
                ));
            }
        },
        Expr::Len(x) => {
            match eval(x, rec)? {
                Value::List(l) => Value::Int(l.len() as i64),
                other => Value::Int(other.as_bytes().map(|b| b.len() as i64).ok_or_else(|| {
                    KiraError::Expression(format!("len() of {}", other.type_name()))
                })?),
            }
        }
        Expr::Min(items) | Expr::Max(items) => {
            let want_min = matches!(e, Expr::Min(_));
            let mut best: Option<Value<'a>> = None;
            for it in items {
                let v = eval(it, rec)?;
                best = Some(match best {
                    None => v,
                    Some(b) => {
                        let take = compare(if want_min { BinOp::Lt } else { BinOp::Gt }, &v, &b)?;
                        if take { v } else { b }
                    }
                });
            }
            best.unwrap_or(Value::Bool(false))
        }
        Expr::CsvMatch(x, set) => {
            let v = eval(x, rec)?;
            let s = v
                .as_bytes()
                .ok_or_else(|| KiraError::Expression("csv_match() requires a string".into()))?;
            Value::Bool(set.iter().any(|i| i == s))
        }
        Expr::WildcardMatch(x, re) | Expr::RegexMatch(x, re) => {
            let v = eval(x, rec)?;
            let s = v
                .as_bytes()
                .ok_or_else(|| KiraError::Expression("match functions require a string".into()))?;
            Value::Bool(re.is_match(s))
        }
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::pairs::record::split_fields;

    fn cols() -> ColumnMap {
        let names: Vec<String> = [
            "readID",
            "chrom1",
            "pos1",
            "chrom2",
            "pos2",
            "strand1",
            "strand2",
            "pair_type",
            "mapq1",
            "mapq2",
            "phase1",
        ]
        .iter()
        .map(|s| s.to_string())
        .collect();
        ColumnMap::from_names(names).unwrap()
    }

    fn run(expr: &str, line: &str) -> Result<bool> {
        let c = cols();
        let f = Filter::compile(expr, &c, &[])?;
        let mut ends = Vec::new();
        split_fields(line.as_bytes(), &mut ends);
        f.eval(&PairRecordRef::new(line.as_bytes(), &ends))
    }

    const L: &str = "r1\tchr1\t100\tchr1\t1500\t+\t-\tUU\t60\t10\t.";
    const T: &str = "r2\tchr1\t100\tchr2\t1500\t+\t-\tUR\t60\t40\t1";

    #[test]
    fn evaluates_examples() {
        assert!(run("pair_type == \"UU\"", L).unwrap());
        assert!(!run("pair_type == \"UU\"", T).unwrap());
        assert!(run("(pair_type == \"UU\") and (abs(pos1-pos2) < 2000)", L).unwrap());
        assert!(!run("(pair_type == \"UU\") and (abs(pos1-pos2) < 1000)", L).unwrap());
        assert!(run("chrom1 == chrom2", L).unwrap());
        assert!(!run("chrom1 == chrom2", T).unwrap());
        assert!(run("mapq1 >= 30 and mapq2 >= 30", T).unwrap());
        assert!(!run("mapq1 >= 30 and mapq2 >= 30", L).unwrap());
        assert!(
            run(
                "(pair_type==\"UU\") or (pair_type==\"UR\") or (pair_type==\"RU\")",
                T
            )
            .unwrap()
        );
        assert!(run("COLS[1]==COLS[3]", L).unwrap());
        assert!(run("(chrom1==chrom2) and (abs(pos1 - pos2) < 1e6)", L).unwrap());
        assert!(!run("(chrom1==\"!\") and (chrom2!=\"!\")", L).unwrap());
        assert!(
            run(
                "regex_match(chrom1, \"chr\\d+\") and regex_match(chrom2, \"chr\\d+\")",
                T
            )
            .unwrap()
        );
        assert!(run("wildcard_match(pair_type, 'U*')", T).unwrap());
        assert!(!run("wildcard_match(pair_type, 'R*')", T).unwrap());
        assert!(run("csv_match(chrom2, 'chr2,chr3')", T).unwrap());
        assert!(run("pair_type in ['UU', 'UR']", T).unwrap());
        assert!(run("chrom2 not in ['chr2']", L).unwrap());
        assert!(run("True", L).unwrap());
        assert!(run("not (pos1 > pos2)", L).unwrap());
        assert!(run("phase1 == '.'", L).unwrap());
        assert!(run("int(phase1) == 1", T).unwrap());
        assert!(run("(pos2 - pos1) // 1000 == 1", L).unwrap());
        assert!(run("(pos2 - pos1) % 1000 == 400", L).unwrap());
        assert!(run("1 < mapq2 < 20", L).unwrap());
        assert!(run("min(mapq1, mapq2) == 10", L).unwrap());
        assert!(run("-pos1 == -100", L).unwrap());
        assert!(run("7 // -2 == -4 and -7 // 2 == -4 and 7 // 2 == 3", L).unwrap());
        assert!(run("7 % -2 == -1 and -7 % 2 == 1 and 7 % 2 == 1", L).unwrap());
    }

    #[test]
    fn errors() {
        assert!(Filter::compile("nope == 1", &cols(), &[]).is_err());
        assert!(run("phase1 > 5", T).is_err());
        assert!(run("int(phase1) == 1", L).is_err());
        assert!(Filter::compile("foo(1)", &cols(), &[]).is_err());
        let f =
            Filter::compile("phase1 > 0", &cols(), &[("phase1".into(), ColumnType::Int)]).unwrap();
        let mut ends = Vec::new();
        split_fields(T.as_bytes(), &mut ends);
        assert!(f.eval(&PairRecordRef::new(T.as_bytes(), &ends)).unwrap());
    }
}
