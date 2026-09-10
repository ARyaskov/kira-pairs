//! Minimal stderr logger driven by `-v` flags.

use log::{Level, LevelFilter, Log, Metadata, Record};
use std::io::Write;
use std::time::Instant;

struct StderrLogger {
    start: Instant,
}

impl Log for StderrLogger {
    fn enabled(&self, metadata: &Metadata) -> bool {
        metadata.level() <= log::max_level()
    }

    fn log(&self, record: &Record) {
        if !self.enabled(record.metadata()) {
            return;
        }
        let elapsed = self.start.elapsed().as_secs_f64();
        let level = match record.level() {
            Level::Error => "error",
            Level::Warn => "warning",
            Level::Info => "info",
            Level::Debug => "debug",
            Level::Trace => "trace",
        };
        let stderr = std::io::stderr();
        let mut lock = stderr.lock();
        // Ignore failures: diagnostics must never abort data processing.
        let _ = writeln!(
            lock,
            "[kira-pairs {elapsed:8.2}s {level}] {}",
            record.args()
        );
    }

    fn flush(&self) {
        let _ = std::io::stderr().flush();
    }
}

/// Install the global logger. `verbosity` is the number of `-v` flags.
pub fn init(verbosity: u8, quiet: bool) {
    let level = if quiet {
        LevelFilter::Error
    } else {
        match verbosity {
            0 => LevelFilter::Warn,
            1 => LevelFilter::Info,
            2 => LevelFilter::Debug,
            _ => LevelFilter::Trace,
        }
    };
    let logger = Box::new(StderrLogger {
        start: Instant::now(),
    });
    // A second initialisation (e.g. in tests) is harmless: keep the first
    // logger but always apply the requested level.
    let _ = log::set_boxed_logger(logger);
    log::set_max_level(level);
}
