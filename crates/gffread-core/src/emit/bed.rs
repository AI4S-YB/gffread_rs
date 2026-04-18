use std::io::{self, Write};

use crate::model::{Annotation, Transcript};

pub fn write_bed<W: Write>(out: &mut W, annotation: &Annotation) -> io::Result<()> {
    for transcript in &annotation.transcripts {
        write_transcript(out, transcript)?;
    }

    Ok(())
}

fn write_transcript<W: Write>(out: &mut W, transcript: &Transcript) -> io::Result<()> {
    let thick_start = transcript
        .cds
        .first()
        .map_or(transcript.start, |segment| segment.start)
        - 1;
    let thick_end = transcript
        .cds
        .last()
        .map_or(transcript.end, |segment| segment.end);
    let item_rgb = format!("{},0,0", transcript.cds_phase().unwrap_or('0'));
    let block_count = transcript.exons.len().max(1);
    let block_sizes = if transcript.exons.is_empty() {
        format!("{},", transcript.end - transcript.start + 1)
    } else {
        transcript
            .exons
            .iter()
            .map(|segment| format!("{},", segment.end - segment.start + 1))
            .collect::<String>()
    };
    let block_starts = if transcript.exons.is_empty() {
        "0,".to_owned()
    } else {
        transcript
            .exons
            .iter()
            .map(|segment| format!("{},", segment.start - transcript.start))
            .collect::<String>()
    };

    write!(
        out,
        "{}\t{}\t{}\t{}\t100\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
        transcript.seqid,
        transcript.start - 1,
        transcript.end,
        transcript.id,
        transcript.strand,
        thick_start,
        thick_end,
        item_rgb,
        block_count,
        block_sizes,
        block_starts
    )?;

    let mut wrote_attr = false;
    if transcript.has_cds() {
        write!(
            out,
            "\tCDS={}:{}",
            transcript.cds.first().expect("checked has_cds").start - 1,
            transcript.cds.last().expect("checked has_cds").end
        )?;
        wrote_attr = true;
    }
    if let Some(phase) = transcript.cds_phase() {
        if wrote_attr {
            write!(out, ";")?;
        } else {
            write!(out, "\t")?;
        }
        write!(out, "CDSphase={phase}")?;
        wrote_attr = true;
    }
    if let Some(gene_id) = &transcript.gene_id {
        if wrote_attr {
            write!(out, ";")?;
        } else {
            write!(out, "\t")?;
        }
        write!(out, "geneID={gene_id}")?;
        wrote_attr = true;
    }
    if let Some(gene_name) = &transcript.gene_name {
        if wrote_attr {
            write!(out, ";")?;
        } else {
            write!(out, "\t")?;
        }
        write!(out, "gene_name={gene_name}")?;
    }

    writeln!(out)
}
