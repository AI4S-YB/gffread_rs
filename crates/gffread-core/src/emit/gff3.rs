use std::io::{self, Write};

use crate::model::{Annotation, Transcript};

pub fn write_gff3<W: Write>(
    out: &mut W,
    annotation: &Annotation,
    version: &str,
    command_line: &str,
) -> io::Result<()> {
    writeln!(out, "##gff-version 3")?;
    writeln!(out, "# gffread v{version}")?;
    writeln!(out, "# {command_line}")?;

    for transcript in &annotation.transcripts {
        write_transcript(out, transcript)?;
    }

    Ok(())
}

fn write_transcript<W: Write>(out: &mut W, transcript: &Transcript) -> io::Result<()> {
    write!(
        out,
        "{}\t{}\t{}\t{}\t{}\t.\t{}\t.\tID={}",
        transcript.seqid,
        transcript.source,
        transcript.feature,
        transcript.start,
        transcript.end,
        transcript.strand,
        transcript.id
    )?;

    if let Some(gene_id) = &transcript.gene_id {
        write!(out, ";geneID={gene_id}")?;
    }

    if let Some(gene_name) = &transcript.gene_name {
        write!(out, ";gene_name={gene_name}")?;
    }

    writeln!(out)?;

    for exon in &transcript.exons {
        writeln!(
            out,
            "{}\t{}\texon\t{}\t{}\t.\t{}\t.\tParent={}",
            transcript.seqid,
            transcript.source,
            exon.start,
            exon.end,
            transcript.strand,
            transcript.id
        )?;
    }

    for cds in &transcript.cds {
        writeln!(
            out,
            "{}\t{}\tCDS\t{}\t{}\t.\t{}\t{}\tParent={}",
            transcript.seqid,
            transcript.source,
            cds.start,
            cds.end,
            transcript.strand,
            cds.phase,
            transcript.id
        )?;
    }

    Ok(())
}
