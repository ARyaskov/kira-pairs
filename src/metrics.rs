//! Optional runtime metrics (`--metrics`) and periodic progress
//! (`--progress`) reporting. Both are designed to add negligible overhead
//! when disabled: counters are plain integers owned by the pipeline and the
//! progress check is a cheap modulo test.

use std::collections::BTreeMap;
use std::fmt::Write as _;
use std::io::Write;
use std::time::{Duration, Instant};

/// Job-level counters printed at exit when `--metrics` is set.
#[derive(Debug, Default, Clone)]
pub struct Metrics {
    values: BTreeMap<&'static str, MetricValue>,
    start: Option<Instant>,
}

/// A metric value.
#[derive(Debug, Clone)]
pub enum MetricValue {
    /// Counter.
    Int(u64),
    /// Duration or rate.
    Float(f64),
}

impl Metrics {
    /// Start the wall clock.
    pub fn start() -> Self {
        Self {
            values: BTreeMap::new(),
            start: Some(Instant::now()),
        }
    }

    /// Set an integer metric.
    pub fn set(&mut self, key: &'static str, v: u64) {
        self.values.insert(key, MetricValue::Int(v));
    }

    /// Set a float metric (seconds, rates).
    pub fn set_f(&mut self, key: &'static str, v: f64) {
        self.values.insert(key, MetricValue::Float(v));
    }

    /// Add to an integer metric.
    pub fn add(&mut self, key: &'static str, v: u64) {
        match self.values.get_mut(key) {
            Some(MetricValue::Int(x)) => *x += v,
            _ => {
                self.values.insert(key, MetricValue::Int(v));
            }
        }
    }

    /// Seconds since [`Metrics::start`].
    pub fn elapsed(&self) -> f64 {
        self.start.map(|s| s.elapsed().as_secs_f64()).unwrap_or(0.0)
    }

    /// Finalise derived metrics (elapsed, rate, peak RSS).
    pub fn finalize(&mut self) {
        let secs = self.elapsed();
        self.set_f("wall_seconds", secs);
        if let Some(MetricValue::Int(n)) = self.values.get("records_read")
            && secs > 0.0
        {
            let rate = *n as f64 / secs;
            self.set_f("records_per_second", rate);
        }
        if let Some(rss) = peak_rss_bytes() {
            self.set("peak_rss_bytes", rss);
        }
    }

    /// Render as `key\tvalue` lines.
    pub fn render(&self) -> String {
        let mut out = String::new();
        for (k, v) in &self.values {
            match v {
                MetricValue::Int(i) => {
                    let _ = writeln!(out, "{k}\t{i}");
                }
                MetricValue::Float(f) => {
                    let _ = writeln!(out, "{k}\t{f:.3}");
                }
            }
        }
        out
    }

    /// Print to stderr.
    pub fn print(&self) {
        let stderr = std::io::stderr();
        let mut l = stderr.lock();
        let _ = writeln!(l, "# kira-pairs metrics");
        let _ = l.write_all(self.render().as_bytes());
    }
}

/// Peak resident set size in bytes (Linux: `VmHWM` from `/proc/self/status`).
pub fn peak_rss_bytes() -> Option<u64> {
    let status = std::fs::read_to_string("/proc/self/status").ok()?;
    for line in status.lines() {
        if let Some(rest) = line.strip_prefix("VmHWM:") {
            let kb: u64 = rest.trim().trim_end_matches("kB").trim().parse().ok()?;
            return Some(kb * 1024);
        }
    }
    None
}

/// Current resident set size in bytes.
pub fn current_rss_bytes() -> Option<u64> {
    let status = std::fs::read_to_string("/proc/self/status").ok()?;
    for line in status.lines() {
        if let Some(rest) = line.strip_prefix("VmRSS:") {
            let kb: u64 = rest.trim().trim_end_matches("kB").trim().parse().ok()?;
            return Some(kb * 1024);
        }
    }
    None
}

/// Low-overhead periodic progress reporter.
pub struct Progress {
    enabled: bool,
    start: Instant,
    last: Instant,
    interval: Duration,
    last_records: u64,
    tick: u64,
}

impl Progress {
    /// Reporter printing at most every `interval`.
    pub fn new(enabled: bool, interval: Duration) -> Self {
        let now = Instant::now();
        Self {
            enabled,
            start: now,
            last: now,
            interval,
            last_records: 0,
            tick: 0,
        }
    }

    /// Whether progress output is enabled.
    pub fn enabled(&self) -> bool {
        self.enabled
    }

    /// Report progress; cheap when disabled or between intervals.
    #[inline]
    pub fn tick(&mut self, records: u64, bytes_read: u64) {
        if !self.enabled {
            return;
        }
        self.tick += 1;
        if self.tick & 0xffff != 0 {
            return;
        }
        let now = Instant::now();
        if now.duration_since(self.last) < self.interval {
            return;
        }
        self.report(records, bytes_read, now);
    }

    /// Force a report (e.g. at the end).
    pub fn finish(&mut self, records: u64, bytes_read: u64) {
        if self.enabled {
            self.report(records, bytes_read, Instant::now());
        }
    }

    fn report(&mut self, records: u64, bytes_read: u64, now: Instant) {
        let dt = now.duration_since(self.last).as_secs_f64().max(1e-9);
        let rate = (records.saturating_sub(self.last_records)) as f64 / dt;
        let elapsed = now.duration_since(self.start).as_secs_f64();
        let rss = current_rss_bytes().unwrap_or(0);
        let stderr = std::io::stderr();
        let mut l = stderr.lock();
        let _ = writeln!(
            l,
            "[kira-pairs progress {elapsed:8.1}s] {:.1} M pairs  {}  {:.2} M pairs/s  {} RSS",
            records as f64 / 1e6,
            crate::memory::format_size(bytes_read),
            rate / 1e6,
            crate::memory::format_size(rss)
        );
        self.last = now;
        self.last_records = records;
    }
}
