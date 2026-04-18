use std::io::{self, Write};

use crate::model::{Annotation, Transcript};

pub fn write_tlf<W: Write>(
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
    write!(
        out,
        "{}\t{}\t{}\t{}\t{}\t.\t{}\t.\tID={};exonCount={};exons={}",
        transcript.seqid,
        source,
        transcript.feature,
        transcript.start,
        transcript.end,
        transcript.strand,
        transcript.id,
        transcript.exons.len(),
        transcript.exon_list()
    )?;

    if transcript.has_cds() {
        let cds_start = transcript.cds.first().expect("checked has_cds").start;
        let cds_end = transcript.cds.last().expect("checked has_cds").end;
        write!(out, ";CDS={cds_start}:{cds_end}")?;
    }

    if let Some(phase) = transcript.cds_phase() {
        write!(out, ";CDSphase={phase}")?;
    }

    if let Some(gene_id) = &transcript.gene_id {
        write!(out, ";geneID={gene_id}")?;
    }

    if let Some(gene_name) = &transcript.gene_name {
        write!(out, ";gene_name={gene_name}")?;
    }

    writeln!(out)
}
