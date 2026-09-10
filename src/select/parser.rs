//! Tokenizer and recursive-descent parser producing an AST.

use crate::error::{KiraError, Result};

/// Binary operators.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum BinOp {
    /// `==`
    Eq,
    /// `!=`
    Ne,
    /// `<`
    Lt,
    /// `<=`
    Le,
    /// `>`
    Gt,
    /// `>=`
    Ge,
    /// `in`
    In,
    /// `not in`
    NotIn,
    /// `+`
    Add,
    /// `-`
    Sub,
    /// `*`
    Mul,
    /// `/`
    Div,
    /// `//`
    FloorDiv,
    /// `%`
    Mod,
    /// `**`
    Pow,
}

/// Abstract syntax tree.
#[derive(Debug, Clone, PartialEq)]
pub enum Ast {
    /// Integer literal.
    Int(i64),
    /// Float literal.
    Float(f64),
    /// String literal.
    Str(String),
    /// Boolean literal.
    Bool(bool),
    /// Identifier (column name or constant).
    Name(String),
    /// `COLS[i]`
    ColIndex(usize),
    /// List literal.
    List(Vec<Ast>),
    /// Unary negation.
    Neg(Box<Ast>),
    /// Logical not.
    Not(Box<Ast>),
    /// Logical and.
    And(Box<Ast>, Box<Ast>),
    /// Logical or.
    Or(Box<Ast>, Box<Ast>),
    /// Binary operator.
    Bin(BinOp, Box<Ast>, Box<Ast>),
    /// Function call.
    Call(String, Vec<Ast>),
}

#[derive(Debug, Clone, PartialEq)]
enum Tok {
    Int(i64),
    Float(f64),
    Str(String),
    Name(String),
    Op(&'static str),
    LParen,
    RParen,
    LBracket,
    RBracket,
    Comma,
    End,
}

fn tokenize(src: &str) -> Result<Vec<Tok>> {
    let b = src.as_bytes();
    let mut i = 0;
    let mut out = Vec::new();
    let err = |msg: String| KiraError::Expression(msg);
    while i < b.len() {
        let c = b[i];
        if c.is_ascii_whitespace() {
            i += 1;
            continue;
        }
        if c.is_ascii_digit() || (c == b'.' && i + 1 < b.len() && b[i + 1].is_ascii_digit()) {
            let start = i;
            let mut is_float = false;
            while i < b.len() && (b[i].is_ascii_digit() || b[i] == b'.' || b[i] == b'_') {
                if b[i] == b'.' {
                    is_float = true;
                }
                i += 1;
            }
            if i < b.len() && (b[i] == b'e' || b[i] == b'E') {
                let mut j = i + 1;
                if j < b.len() && (b[j] == b'+' || b[j] == b'-') {
                    j += 1;
                }
                if j < b.len() && b[j].is_ascii_digit() {
                    is_float = true;
                    i = j;
                    while i < b.len() && b[i].is_ascii_digit() {
                        i += 1;
                    }
                }
            }
            let text: String = src[start..i].chars().filter(|c| *c != '_').collect();
            if is_float {
                out.push(Tok::Float(
                    text.parse()
                        .map_err(|_| err(format!("bad number {text:?}")))?,
                ));
            } else {
                out.push(Tok::Int(
                    text.parse()
                        .map_err(|_| err(format!("bad integer {text:?}")))?,
                ));
            }
            continue;
        }
        if c == b'"' || c == b'\'' {
            let quote = c;
            i += 1;
            let mut s = String::new();
            loop {
                if i >= b.len() {
                    return Err(err("unterminated string literal".into()));
                }
                if b[i] == b'\\' && i + 1 < b.len() {
                    let e = b[i + 1];
                    match e {
                        b'n' => s.push('\n'),
                        b't' => s.push('\t'),
                        b'\\' => s.push('\\'),
                        b'\'' => s.push('\''),
                        b'"' => s.push('"'),
                        // Unknown escapes are kept verbatim, as Python does
                        // (needed for regex patterns such as "chr\d+").
                        other => {
                            s.push('\\');
                            s.push(other as char);
                        }
                    }
                    i += 2;
                    continue;
                }
                if b[i] == quote {
                    i += 1;
                    break;
                }
                let ch = src[i..].chars().next().unwrap_or('?');
                s.push(ch);
                i += ch.len_utf8();
            }
            out.push(Tok::Str(s));
            continue;
        }
        if c.is_ascii_alphabetic() || c == b'_' {
            let start = i;
            while i < b.len() && (b[i].is_ascii_alphanumeric() || b[i] == b'_') {
                i += 1;
            }
            out.push(Tok::Name(src[start..i].to_string()));
            continue;
        }
        let two = if i + 1 < b.len() { &src[i..i + 2] } else { "" };
        let op2 = match two {
            "==" => Some("=="),
            "!=" => Some("!="),
            "<=" => Some("<="),
            ">=" => Some(">="),
            "//" => Some("//"),
            "**" => Some("**"),
            "&&" => Some("and"),
            "||" => Some("or"),
            _ => None,
        };
        if let Some(op) = op2 {
            out.push(Tok::Op(op));
            i += 2;
            continue;
        }
        let one = match c {
            b'<' => Some("<"),
            b'>' => Some(">"),
            b'+' => Some("+"),
            b'-' => Some("-"),
            b'*' => Some("*"),
            b'/' => Some("/"),
            b'%' => Some("%"),
            b'!' => Some("not"),
            _ => None,
        };
        if let Some(op) = one {
            out.push(Tok::Op(op));
            i += 1;
            continue;
        }
        match c {
            b'(' => out.push(Tok::LParen),
            b')' => out.push(Tok::RParen),
            b'[' => out.push(Tok::LBracket),
            b']' => out.push(Tok::RBracket),
            b',' => out.push(Tok::Comma),
            b'=' => return Err(err("single '=' is not an operator; use '=='".into())),
            _ => return Err(err(format!("unexpected character {:?}", c as char))),
        }
        i += 1;
    }
    out.push(Tok::End);
    Ok(out)
}

struct Parser {
    toks: Vec<Tok>,
    pos: usize,
}

impl Parser {
    fn peek(&self) -> &Tok {
        &self.toks[self.pos]
    }

    fn next(&mut self) -> Tok {
        let t = self.toks[self.pos].clone();
        if self.pos + 1 < self.toks.len() {
            self.pos += 1;
        }
        t
    }

    fn is_name(&self, n: &str) -> bool {
        matches!(self.peek(), Tok::Name(x) if x == n)
    }

    fn is_op(&self, o: &str) -> bool {
        matches!(self.peek(), Tok::Op(x) if *x == o)
    }

    fn expect(&mut self, t: Tok) -> Result<()> {
        if *self.peek() == t {
            self.next();
            Ok(())
        } else {
            Err(KiraError::Expression(format!(
                "expected {t:?}, found {:?}",
                self.peek()
            )))
        }
    }

    fn or_expr(&mut self) -> Result<Ast> {
        let mut left = self.and_expr()?;
        while self.is_name("or") || self.is_op("or") {
            self.next();
            let right = self.and_expr()?;
            left = Ast::Or(Box::new(left), Box::new(right));
        }
        Ok(left)
    }

    fn and_expr(&mut self) -> Result<Ast> {
        let mut left = self.not_expr()?;
        while self.is_name("and") || self.is_op("and") {
            self.next();
            let right = self.not_expr()?;
            left = Ast::And(Box::new(left), Box::new(right));
        }
        Ok(left)
    }

    fn not_expr(&mut self) -> Result<Ast> {
        if self.is_name("not") || self.is_op("not") {
            self.next();
            let e = self.not_expr()?;
            return Ok(Ast::Not(Box::new(e)));
        }
        self.comparison()
    }

    fn comparison(&mut self) -> Result<Ast> {
        let first = self.arith()?;
        let mut terms = vec![first];
        let mut ops = Vec::new();
        loop {
            let op = match self.peek() {
                Tok::Op("==") => BinOp::Eq,
                Tok::Op("!=") => BinOp::Ne,
                Tok::Op("<") => BinOp::Lt,
                Tok::Op("<=") => BinOp::Le,
                Tok::Op(">") => BinOp::Gt,
                Tok::Op(">=") => BinOp::Ge,
                Tok::Name(n) if n == "in" => BinOp::In,
                Tok::Name(n) if n == "not" => {
                    if matches!(self.toks.get(self.pos + 1), Some(Tok::Name(x)) if x == "in") {
                        self.next();
                        BinOp::NotIn
                    } else {
                        break;
                    }
                }
                _ => break,
            };
            self.next();
            ops.push(op);
            terms.push(self.arith()?);
        }
        if ops.is_empty() {
            return Ok(terms.pop().unwrap_or(Ast::Bool(true)));
        }
        // Python chained comparisons: a < b < c == (a < b) and (b < c).
        let mut result: Option<Ast> = None;
        for (i, op) in ops.into_iter().enumerate() {
            let cmp = Ast::Bin(
                op,
                Box::new(terms[i].clone()),
                Box::new(terms[i + 1].clone()),
            );
            result = Some(match result {
                None => cmp,
                Some(r) => Ast::And(Box::new(r), Box::new(cmp)),
            });
        }
        Ok(result.unwrap_or(Ast::Bool(true)))
    }

    fn arith(&mut self) -> Result<Ast> {
        let mut left = self.term()?;
        loop {
            let op = match self.peek() {
                Tok::Op("+") => BinOp::Add,
                Tok::Op("-") => BinOp::Sub,
                _ => break,
            };
            self.next();
            let right = self.term()?;
            left = Ast::Bin(op, Box::new(left), Box::new(right));
        }
        Ok(left)
    }

    fn term(&mut self) -> Result<Ast> {
        let mut left = self.unary()?;
        loop {
            let op = match self.peek() {
                Tok::Op("*") => BinOp::Mul,
                Tok::Op("/") => BinOp::Div,
                Tok::Op("//") => BinOp::FloorDiv,
                Tok::Op("%") => BinOp::Mod,
                _ => break,
            };
            self.next();
            let right = self.unary()?;
            left = Ast::Bin(op, Box::new(left), Box::new(right));
        }
        Ok(left)
    }

    fn unary(&mut self) -> Result<Ast> {
        if self.is_op("-") {
            self.next();
            let e = self.unary()?;
            return Ok(Ast::Neg(Box::new(e)));
        }
        if self.is_op("+") {
            self.next();
            return self.unary();
        }
        self.power()
    }

    fn power(&mut self) -> Result<Ast> {
        let base = self.atom()?;
        if self.is_op("**") {
            self.next();
            let exp = self.unary()?;
            return Ok(Ast::Bin(BinOp::Pow, Box::new(base), Box::new(exp)));
        }
        Ok(base)
    }

    fn atom(&mut self) -> Result<Ast> {
        match self.next() {
            Tok::Int(i) => Ok(Ast::Int(i)),
            Tok::Float(f) => Ok(Ast::Float(f)),
            Tok::Str(s) => Ok(Ast::Str(s)),
            Tok::LParen => {
                let e = self.or_expr()?;
                self.expect(Tok::RParen)?;
                Ok(e)
            }
            Tok::LBracket => {
                let mut items = Vec::new();
                if *self.peek() != Tok::RBracket {
                    loop {
                        items.push(self.or_expr()?);
                        if *self.peek() == Tok::Comma {
                            self.next();
                            if *self.peek() == Tok::RBracket {
                                break;
                            }
                            continue;
                        }
                        break;
                    }
                }
                self.expect(Tok::RBracket)?;
                Ok(Ast::List(items))
            }
            Tok::Name(n) => {
                match n.as_str() {
                    "True" => return Ok(Ast::Bool(true)),
                    "False" => return Ok(Ast::Bool(false)),
                    _ => {}
                }
                if *self.peek() == Tok::LParen {
                    self.next();
                    let mut args = Vec::new();
                    if *self.peek() != Tok::RParen {
                        loop {
                            args.push(self.or_expr()?);
                            if *self.peek() == Tok::Comma {
                                self.next();
                                continue;
                            }
                            break;
                        }
                    }
                    self.expect(Tok::RParen)?;
                    return Ok(Ast::Call(n, args));
                }
                if n == "COLS" && *self.peek() == Tok::LBracket {
                    self.next();
                    let idx = match self.next() {
                        Tok::Int(i) if i >= 0 => i as usize,
                        other => {
                            return Err(KiraError::Expression(format!(
                                "COLS index must be a non-negative integer, found {other:?}"
                            )));
                        }
                    };
                    self.expect(Tok::RBracket)?;
                    return Ok(Ast::ColIndex(idx));
                }
                Ok(Ast::Name(n))
            }
            other => Err(KiraError::Expression(format!("unexpected token {other:?}"))),
        }
    }
}

/// Parse an expression into an AST.
pub fn parse(src: &str) -> Result<Ast> {
    let toks = tokenize(src)?;
    let mut p = Parser { toks, pos: 0 };
    let ast = p.or_expr()?;
    if *p.peek() != Tok::End {
        return Err(KiraError::Expression(format!(
            "unexpected trailing token {:?}",
            p.peek()
        )));
    }
    Ok(ast)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parses_common_expressions() {
        let a = parse("pair_type == \"UU\"").unwrap();
        assert_eq!(
            a,
            Ast::Bin(
                BinOp::Eq,
                Box::new(Ast::Name("pair_type".into())),
                Box::new(Ast::Str("UU".into()))
            )
        );
        let a = parse("(pair_type == \"UU\") and (abs(pos1-pos2) < 1000)").unwrap();
        assert!(matches!(a, Ast::And(_, _)));
        let a = parse("mapq1 >= 30 and mapq2 >= 30 or not chrom1 == chrom2").unwrap();
        assert!(matches!(a, Ast::Or(_, _)));
        let a = parse("pair_type in ['UU', 'UR', 'RU']").unwrap();
        assert!(matches!(a, Ast::Bin(BinOp::In, _, _)));
        let a = parse("chrom1 not in ['chrM']").unwrap();
        assert!(matches!(a, Ast::Bin(BinOp::NotIn, _, _)));
        let a = parse("COLS[1] == COLS[3]").unwrap();
        assert!(matches!(a, Ast::Bin(BinOp::Eq, _, _)));
        let a = parse("1 < pos1 < 1e6").unwrap();
        assert!(matches!(a, Ast::And(_, _)));
        assert!(parse("pos1 = 5").is_err());
        assert!(parse("(pos1").is_err());
        assert!(parse("pos1 <").is_err());
    }
}
