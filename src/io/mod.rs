//! I/O building blocks: large block reads, compression auto-detection,
//! a parallel BGZF codec and temporary-file management.

pub mod bgzf;
pub mod buffered;
pub mod compression;
pub mod temp;

pub use compression::{
    Compression, FinishWrite, InputSource, open_input, open_output, output_compression_for_path,
};
