use std::io::{self, Write};

use crate::model::{Annotation, Segment, Transcript};

pub fn write_gtf<W: Write>(
    out: &mut W,
    annotation: &Annotation,
    track_label: Option<&str>,
) -> io::Result<()> {
    for transcript in &annotation.transcripts {
        write_transcript(out, transcript, track_label)?;
    }

    Ok(())
}

fn write_transcript<W: Write>(
    out: &mut W,
    transcript: &Transcript,
    track_label: Option<&str>,
) -> io::Result<()> {
    let source = track_label.unwrap_or(&transcript.source);
    writeln!(
        out,
        "{}\t{}\ttranscript\t{}\t{}\t.\t{}\t.\t{}",
        transcript.seqid,
        source,
        transcript.start,
        transcript.end,
        transcript.strand,
        attrs_for_transcript(transcript, false)
    )?;

    for exon in &transcript.exons {
        write_segment(out, transcript, exon, "exon", source)?;
    }

    for cds in &transcript.cds {
        write_segment(out, transcript, cds, "CDS", source)?;
    }

    Ok(())
}

fn write_segment<W: Write>(
    out: &mut W,
    transcript: &Transcript,
    segment: &Segment,
    feature: &str,
    source: &str,
) -> io::Result<()> {
    let phase = if feature == "CDS" {
        segment.phase.as_str()
    } else {
        "."
    };
    let score = if feature == "exon" {
        segment.score.as_str()
    } else {
        "."
    };
    writeln!(
        out,
        "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
        transcript.seqid,
        source,
        feature,
        segment.start,
        segment.end,
        score,
        transcript.strand,
        phase,
        attrs_for_transcript(transcript, true)
    )
}

fn attrs_for_transcript(transcript: &Transcript, terminal_semicolon: bool) -> String {
    let gene_id = transcript.gene_id.as_deref().unwrap_or("");
    let mut attrs = format!(
        "transcript_id \"{}\"; gene_id \"{}\"",
        transcript.id, gene_id
    );

    if let Some(gene_name) = transcript
        .gene_name
        .as_deref()
        .filter(|gene_name| !gene_name.is_empty())
    {
        attrs.push_str(&format!("; gene_name \"{gene_name}\""));
    }
    if terminal_semicolon {
        attrs.push(';');
    }
    attrs
}
