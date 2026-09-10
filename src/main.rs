//! Thin CLI wrapper around the `kira_pairs` library.

use std::process::ExitCode;

use clap::Parser;

fn main() -> ExitCode {
    // Restore default SIGPIPE handling so that a downstream consumer that
    // exits early (e.g. `head`) terminates us quietly, like other Unix
    // tools. This is the only unsafe call in the binary.
    #[cfg(unix)]
    {
        // SAFETY: `signal` with SIG_DFL only resets the disposition of SIGPIPE;
        // it has no memory-safety preconditions and runs before any threads exist.
        unsafe {
            libc::signal(libc::SIGPIPE, libc::SIG_DFL);
        }
    }
    let cli = kira_pairs::cli::Cli::parse();
    kira_pairs::logging::init(cli.verbose, cli.quiet);
    match kira_pairs::cli::run(cli) {
        Ok(()) => ExitCode::SUCCESS,
        Err(e) => {
            if e.is_broken_pipe() {
                // Downstream closed the pipe: exit quietly.
                return ExitCode::from(141);
            }
            eprintln!("kira-pairs: error: {e}");
            ExitCode::from(e.exit_code() as u8)
        }
    }
}
