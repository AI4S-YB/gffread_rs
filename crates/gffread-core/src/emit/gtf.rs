use std::io::{self, Write};

use crate::model::{Annotation, Segment, Transcript};

pub fn write_gtf<W: Write>(out: &mut W, annotation: &Annotation) -> io::Result<()> {
    for transcript in &annotation.transcripts {
        write_transcript(out, transcript)?;
    }

    Ok(())
}

fn write_transcript<W: Write>(out: &mut W, transcript: &Transcript) -> io::Result<()> {
    writeln!(
        out,
        "{}\t{}\ttranscript\t{}\t{}\t.\t{}\t.\t{}",
        transcript.seqid,
        transcript.source,
        transcript.start,
        transcript.end,
        transcript.strand,
        attrs_for_transcript(transcript, false)
    )?;

    for exon in &transcript.exons {
        write_segment(out, transcript, exon, "exon")?;
    }

    for cds in &transcript.cds {
        write_segment(out, transcript, cds, "CDS")?;
    }

    Ok(())
}

fn write_segment<W: Write>(
    out: &mut W,
    transcript: &Transcript,
    segment: &Segment,
    feature: &str,
) -> io::Result<()> {
    writeln!(
        out,
        "{}\t{}\t{}\t{}\t{}\t.\t{}\t{}\t{}",
        transcript.seqid,
        transcript.source,
        feature,
        segment.start,
        segment.end,
        transcript.strand,
        segment.phase,
        attrs_for_transcript(transcript, true)
    )
}

fn attrs_for_transcript(transcript: &Transcript, terminal_semicolon: bool) -> String {
    let gene_id = transcript.gene_id.as_deref().unwrap_or("");
    let gene_name = transcript.gene_name.as_deref().unwrap_or("");

    if terminal_semicolon {
        format!(
            "transcript_id \"{}\"; gene_id \"{}\"; gene_name \"{}\";",
            transcript.id, gene_id, gene_name
        )
    } else {
        format!(
            "transcript_id \"{}\"; gene_id \"{}\"; gene_name \"{}\"",
            transcript.id, gene_id, gene_name
        )
    }
}
