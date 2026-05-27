use std::num::NonZero;
use std::path::PathBuf;

use clap::Parser;
use rsomics_common::{CommonFlags, Result, RsomicsError, Tool, ToolMeta};
use rsomics_help::{Example, FlagSpec, HelpSpec, Origin, Section};

use rsomics_bam_targetcut::{TargetcutOpts, targetcut};

pub const META: ToolMeta = ToolMeta {
    name: env!("CARGO_PKG_NAME"),
    version: env!("CARGO_PKG_VERSION"),
};

#[derive(Parser, Debug)]
#[command(
    name = "rsomics-bam-targetcut",
    version,
    about,
    long_about = None,
    disable_help_flag = true
)]
pub struct Cli {
    /// Input BAM file (must be sorted by coordinate).
    pub input: PathBuf,

    /// Output SAM file (default stdout).
    #[arg(short = 'o', long = "output", default_value = "-")]
    output: String,

    /// Minimum base quality.
    #[arg(short = 'Q', long = "min-baseq", default_value_t = 13)]
    min_base_q: u8,

    /// Penalty for entering a target interval (0→1 HMM transition penalty).
    #[arg(short = 'i', long = "in-penalty", default_value_t = 14000)]
    in_penalty: i32,

    /// Emission score for on-target state at zero coverage.
    #[arg(long = "em0")]
    em0: Option<i32>,

    /// Emission score for on-target state at low coverage (depth byte = 0).
    #[arg(long = "em1")]
    em1: Option<i32>,

    /// Emission score for on-target state at high coverage.
    #[arg(long = "em2")]
    em2: Option<i32>,

    /// Reference FASTA (accepted for compatibility; BAQ realignment is not implemented).
    #[arg(short = 'f', long = "reference")]
    reference: Option<PathBuf>,

    #[command(flatten)]
    pub common: CommonFlags,
}

impl Cli {
    pub fn execute(self) -> Result<()> {
        let opts = TargetcutOpts {
            min_base_q: self.min_base_q,
            in_penalty: self.in_penalty,
            em0: self.em0,
            em1: self.em1,
            em2: self.em2,
            reference: self.reference,
        };

        let output_path = (self.output != "-").then(|| PathBuf::from(&self.output));
        let workers =
            NonZero::new(self.common.threads.unwrap_or(1)).unwrap_or(NonZero::<usize>::MIN);
        let stats = targetcut(&self.input, output_path.as_deref(), &opts, workers)?;

        if self.common.json {
            eprintln!(
                "{}",
                serde_json::to_string(&stats)
                    .map_err(|e| RsomicsError::InvalidInput(format!("JSON: {e}")))?
            );
        }

        Ok(())
    }
}

impl Tool for Cli {
    fn meta() -> ToolMeta {
        META
    }

    fn common(&self) -> &CommonFlags {
        &self.common
    }

    fn execute(self) -> Result<()> {
        self.execute()
    }
}

pub static HELP: HelpSpec = HelpSpec {
    name: META.name,
    version: META.version,
    tagline: "Identify amplicon target intervals from pileup depth.",
    origin: Some(Origin {
        upstream: "samtools targetcut",
        upstream_license: "MIT",
        our_license: "MIT OR Apache-2.0",
        paper_doi: None,
    }),
    usage_lines: &["<in.bam> [-Q minQ] [-i inPenalty] [--em0 N] [--em1 N] [--em2 N]"],
    sections: &[Section {
        title: "OPTIONS",
        flags: &[
            FlagSpec {
                short: Some('Q'),
                long: "min-baseq",
                aliases: &[],
                value: Some("INT"),
                type_hint: None,
                required: false,
                default: Some("13"),
                description: "Minimum base quality to include in pileup.",
                why_default: None,
            },
            FlagSpec {
                short: Some('i'),
                long: "in-penalty",
                aliases: &[],
                value: Some("INT"),
                type_hint: None,
                required: false,
                default: Some("14000"),
                description: "HMM cost for transitioning into a target interval.",
                why_default: None,
            },
            FlagSpec {
                short: None,
                long: "em0",
                aliases: &[],
                value: Some("INT"),
                type_hint: None,
                required: false,
                default: Some("-4"),
                description: "On-target emission score for zero-coverage positions.",
                why_default: None,
            },
            FlagSpec {
                short: None,
                long: "em1",
                aliases: &[],
                value: Some("INT"),
                type_hint: None,
                required: false,
                default: Some("1"),
                description: "On-target emission score for low-coverage positions.",
                why_default: None,
            },
            FlagSpec {
                short: None,
                long: "em2",
                aliases: &[],
                value: Some("INT"),
                type_hint: None,
                required: false,
                default: Some("6"),
                description: "On-target emission score for high-coverage positions.",
                why_default: None,
            },
            FlagSpec {
                short: Some('f'),
                long: "reference",
                aliases: &[],
                value: Some("FILE"),
                type_hint: None,
                required: false,
                default: None,
                description: "Reference FASTA (accepted for CLI compatibility; BAQ not implemented).",
                why_default: None,
            },
        ],
    }],
    examples: &[
        Example {
            description: "Find target intervals in an amplicon BAM",
            command: "rsomics-bam-targetcut amplicon.bam",
        },
        Example {
            description: "Raise coverage threshold for calling targets",
            command: "rsomics-bam-targetcut -Q 20 amplicon.bam",
        },
    ],
    json_result_schema_doc: None,
};

#[cfg(test)]
mod tests {
    use super::*;
    use clap::CommandFactory;

    #[test]
    fn cli_debug_assert() {
        Cli::command().debug_assert();
    }
}
