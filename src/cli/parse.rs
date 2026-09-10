//! `kira-pairs parse`.

use std::io::Write;
use std::path::PathBuf;

use clap::Args;

use crate::chroms::ChromSizes;
use crate::cli::common::{Context, Resources, append_pg, new_dict, stdio_path};
use crate::cli::{ResourceOpts, WalksPolicyArg};
use crate::error::Result;
use crate::metrics::{Metrics, Progress};
use crate::parse::hic::WalksPolicy;
use crate::parse::{HicParser, ParseConfig, drive, open_alignments};
use crate::stats::{DistBins, StatsAccumulator, StatsFormat};

/// Parse-related options shared with `process`.
#[derive(Debug, Clone, Args)]
pub struct ParseOpts {
    /// Chromosome order (chrom.sizes-like file). Scaffolds not listed
    /// follow in lexicographic order.
    #[arg(short = 'c', long, value_name = "FILE")]
    pub chroms_path: PathBuf,
    /// Genome assembly name stored in the header.
    #[arg(long, value_name = "NAME")]
    pub assembly: Option<String>,
    /// Minimal MAPQ to consider an alignment uniquely mapped.
    #[arg(long, default_value_t = 1, value_name = "N")]
    pub min_mapq: u8,
    /// Maximal Hi-C molecule size used to rescue single ligations.
    #[arg(long, default_value_t = 750, value_name = "N")]
    pub max_molecule_size: u64,
    /// Read segments not covered by any alignment and longer than this
    /// are treated as null alignments.
    #[arg(long, default_value_t = 20, value_name = "N")]
    pub max_inter_align_gap: u64,
    /// Policy for unrescuable walks.
    #[arg(long, value_enum, default_value_t = WalksPolicyArg::FiveUnique)]
    pub walks_policy: WalksPolicyArg,
    /// Report the 3' end of alignments instead of the 5' end.
    #[arg(long, default_value = "5", value_parser = ["5", "3"], value_name = "END")]
    pub report_alignment_end: String,
    /// Do not flip pairs into upper-triangular order.
    #[arg(long)]
    pub no_flip: bool,
    /// Replace read IDs with `.`.
    #[arg(long)]
    pub drop_readid: bool,
    /// Remove sequences and qualities from SAM columns.
    #[arg(long)]
    pub drop_seq: bool,
    /// Do not output sam1/sam2 columns.
    #[arg(long)]
    pub drop_sam: bool,
    /// Add walk_pair_index/walk_pair_type columns.
    #[arg(long)]
    pub add_pair_index: bool,
    /// Comma-separated extra columns (mapq, pos5, pos3, cigar, read_len,
    /// matched_bp, algn_ref_span, algn_read_span, dist_to_5, dist_to_3, seq,
    /// read_side, algn_idx, same_side_algn_count, or a two-letter SAM tag).
    #[arg(long, value_name = "LIST", default_value = "")]
    pub add_columns: String,
}

impl ParseOpts {
    /// Convert to a library config.
    pub fn to_config(&self) -> ParseConfig {
        ParseConfig {
            min_mapq: self.min_mapq,
            max_molecule_size: self.max_molecule_size,
            max_inter_align_gap: Some(self.max_inter_align_gap),
            walks_policy: match self.walks_policy {
                WalksPolicyArg::Mask => WalksPolicy::Mask,
                WalksPolicyArg::FiveAny => WalksPolicy::FiveAny,
                WalksPolicyArg::FiveUnique => WalksPolicy::FiveUnique,
                WalksPolicyArg::ThreeAny => WalksPolicy::ThreeAny,
                WalksPolicyArg::ThreeUnique => WalksPolicy::ThreeUnique,
            },
            report_3_end: self.report_alignment_end == "3",
            flip: !self.no_flip,
            drop_readid: self.drop_readid,
            drop_seq: self.drop_seq,
            drop_sam: self.drop_sam,
            add_pair_index: self.add_pair_index,
            add_columns: self
                .add_columns
                .split(',')
                .map(str::trim)
                .filter(|s| !s.is_empty())
                .map(String::from)
                .collect(),
            assembly: self.assembly.clone(),
        }
    }
}

/// Arguments for `parse`.
#[derive(Debug, Args)]
pub struct ParseArgs {
    /// Input SAM or BAM (name-grouped, as produced by bwa mem); `-`/omitted = stdin.
    pub input: Option<PathBuf>,
    /// Output .pairs file; `-` or omitted = stdout.
    #[arg(short, long, value_name = "FILE")]
    pub output: Option<PathBuf>,
    /// Output file for statistics.
    #[arg(long, value_name = "FILE")]
    pub output_stats: Option<PathBuf>,
    /// Parse options.
    #[command(flatten)]
    pub parse: ParseOpts,
    /// Resource options.
    #[command(flatten)]
    pub res: ResourceOpts,
}

/// Run `parse`.
pub fn run(a: ParseArgs, ctx: &Context) -> Result<()> {
    let res = Resources::from_opts(&a.res)?;
    let mut metrics = Metrics::start();
    let chroms = ChromSizes::from_path(&a.parse.chroms_path)?;
    let cfg = a.parse.to_config();
    let (source, sam_header, name) = open_alignments(stdio_path(&a.input), res.io_threads)?;
    let dict = new_dict();
    let mut parser = HicParser::new(cfg, sam_header, &chroms, dict.clone())?;
    append_pg(parser.header_mut(), "parse", ctx)?;
    let mut writer = res.open_writer(stdio_path(&a.output))?;
    writer.write_header(parser.header())?;
    let mut stats = a.output_stats.as_ref().map(|_| {
        StatsAccumulator::new(DistBins::default(), dict.get(crate::chroms::UNMAPPED_CHROM))
    });
    let mut progress = Progress::new(ctx.progress, std::time::Duration::from_secs(5));
    let mut n = 0u64;
    drive(source, &mut parser, stats.as_mut(), &name, |_key, line| {
        n += 1;
        progress.tick(n, 0);
        writer.write_line(line)
    })?;
    let bytes_written = writer.finish()?;
    if let (Some(path), Some(s)) = (&a.output_stats, stats) {
        let text = crate::stats::format::render(&s.snapshot(&dict), StatsFormat::Tsv);
        let mut w = crate::io::compression::open_output(
            stdio_path(&Some(path.clone())),
            res.output_options(),
        )?;
        w.write_all(text.as_bytes())?;
        w.finish()?;
    }
    progress.finish(n, 0);
    if ctx.metrics {
        metrics.set("alignments_read", parser.records_in);
        metrics.set("records_written", parser.pairs_out);
        metrics.set("bytes_written", bytes_written);
        metrics.finalize();
        metrics.print();
    }
    Ok(())
}
