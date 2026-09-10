//! Typed error definitions shared by the library and the CLI.

use std::fmt;
use std::path::PathBuf;

/// Location information attached to data errors when it is known.
#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct Location {
    /// Source file name (`-` for stdin) when known.
    pub file: Option<String>,
    /// 1-based line number within the file when known.
    pub line: Option<u64>,
    /// 1-based column index (TAB separated field) when known.
    pub column: Option<usize>,
    /// Offending value, when known and printable.
    pub value: Option<String>,
}

impl Location {
    /// Location with only a file name.
    pub fn file(file: impl Into<String>) -> Self {
        Self {
            file: Some(file.into()),
            ..Default::default()
        }
    }

    /// Attach a line number.
    #[must_use]
    pub fn at_line(mut self, line: u64) -> Self {
        self.line = Some(line);
        self
    }

    /// Attach a column number.
    #[must_use]
    pub fn at_column(mut self, column: usize) -> Self {
        self.column = Some(column);
        self
    }

    /// Attach the offending value (truncated to keep messages readable).
    #[must_use]
    pub fn with_value(mut self, value: &[u8]) -> Self {
        let mut s = String::from_utf8_lossy(value).into_owned();
        if s.len() > 64 {
            s.truncate(61);
            s.push_str("...");
        }
        self.value = Some(s);
        self
    }
}

impl fmt::Display for Location {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let mut parts = Vec::new();
        if let Some(file) = &self.file {
            parts.push(format!("file {file}"));
        }
        if let Some(line) = self.line {
            parts.push(format!("line {line}"));
        }
        if let Some(col) = self.column {
            parts.push(format!("column {col}"));
        }
        if let Some(v) = &self.value {
            parts.push(format!("value {v:?}"));
        }
        if parts.is_empty() {
            Ok(())
        } else {
            write!(f, " ({})", parts.join(", "))
        }
    }
}

/// All errors produced by kira-pairs.
#[derive(Debug, thiserror::Error)]
pub enum KiraError {
    /// I/O failure with the path that was being accessed.
    #[error("I/O error on {path}: {source}")]
    Io {
        /// Path being accessed.
        path: PathBuf,
        /// Underlying error.
        #[source]
        source: std::io::Error,
    },
    /// I/O failure without an associated path (pipes, threads).
    #[error("I/O error: {0}")]
    IoPlain(#[from] std::io::Error),
    /// Malformed `.pairs` data.
    #[error("malformed pairs data{location}: {message}")]
    Format {
        /// What went wrong.
        message: String,
        /// Where it went wrong.
        location: Location,
    },
    /// Malformed or inconsistent header.
    #[error("invalid pairs header{location}: {message}")]
    Header {
        /// What went wrong.
        message: String,
        /// Where it went wrong.
        location: Location,
    },
    /// Malformed chromosome sizes file.
    #[error("invalid chromosome sizes{location}: {message}")]
    ChromSizes {
        /// What went wrong.
        message: String,
        /// Where it went wrong.
        location: Location,
    },
    /// A required column was not found.
    #[error("missing required column {column:?}{location}")]
    MissingColumn {
        /// Column name.
        column: String,
        /// Where it was required.
        location: Location,
    },
    /// Invalid command-line arguments or option combinations.
    #[error("invalid argument: {0}")]
    InvalidArgument(String),
    /// Compressed input could not be decoded (truncated or corrupt).
    #[error("compression error on {path}: {message}")]
    Compression {
        /// Input path.
        path: PathBuf,
        /// Details.
        message: String,
    },
    /// Problems reading SAM/BAM input.
    #[error("alignment input error{location}: {message}")]
    Alignment {
        /// What went wrong.
        message: String,
        /// Where it went wrong.
        location: Location,
    },
    /// Private run file was corrupt or produced by an incompatible version.
    #[error("temporary run file {path} is corrupt or incompatible: {message}")]
    RunFile {
        /// Run file path.
        path: PathBuf,
        /// Details.
        message: String,
    },
    /// Filter expression errors.
    #[error("invalid filter expression: {0}")]
    Expression(String),
    /// Input violated an ordering precondition (e.g. dedup on unsorted pairs).
    #[error("input is not sorted{location}: {message}")]
    NotSorted {
        /// What went wrong.
        message: String,
        /// Where it went wrong.
        location: Location,
    },
    /// Unsupported feature encountered in the input.
    #[error("unsupported: {0}")]
    Unsupported(String),
    /// Memory budget cannot accommodate the request.
    #[error("memory budget too small: {0}")]
    Memory(String),
}

impl KiraError {
    /// Create an I/O error carrying a path.
    pub fn io(path: impl Into<PathBuf>, source: std::io::Error) -> Self {
        Self::Io {
            path: path.into(),
            source,
        }
    }

    /// Create a data format error.
    pub fn format(message: impl Into<String>, location: Location) -> Self {
        Self::Format {
            message: message.into(),
            location,
        }
    }

    /// Create a header error.
    pub fn header(message: impl Into<String>) -> Self {
        Self::Header {
            message: message.into(),
            location: Location::default(),
        }
    }

    /// Create an invalid-argument error.
    pub fn arg(message: impl Into<String>) -> Self {
        Self::InvalidArgument(message.into())
    }

    /// True when the error is a broken pipe (downstream consumer went away).
    pub fn is_broken_pipe(&self) -> bool {
        match self {
            Self::Io { source, .. } => source.kind() == std::io::ErrorKind::BrokenPipe,
            Self::IoPlain(e) => e.kind() == std::io::ErrorKind::BrokenPipe,
            _ => false,
        }
    }

    /// Process exit code for this error class.
    pub fn exit_code(&self) -> i32 {
        match self {
            Self::InvalidArgument(_) | Self::Expression(_) => 2,
            Self::Io { .. } | Self::IoPlain(_) | Self::Compression { .. } => 3,
            Self::Format { .. }
            | Self::Header { .. }
            | Self::ChromSizes { .. }
            | Self::MissingColumn { .. }
            | Self::Alignment { .. }
            | Self::NotSorted { .. } => 4,
            Self::RunFile { .. } | Self::Memory(_) => 5,
            Self::Unsupported(_) => 6,
        }
    }
}

/// Convenience alias used across the crate.
pub type Result<T> = std::result::Result<T, KiraError>;
