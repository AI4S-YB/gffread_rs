use std::io::{self, Write};

use crate::emit::locus::write_locus;
use crate::model::{Annotation, Locus, Transcript};

pub fn write_gff3<W: Write>(
    out: &mut W,
    annotation: &Annotation,
    loci: &[Locus],
    version: &str,
    command_line: &str,
    attrs_filter: Option<&[String]>,
    keep_all_attrs: bool,
) -> io::Result<()> {
    writeln!(out, "##gff-version 3")?;
    writeln!(out, "# gffread v{version}")?;
    writeln!(out, "# {command_line}")?;

    if loci.is_empty() {
        for transcript in &annotation.transcripts {
            write_transcript(out, transcript, attrs_filter, keep_all_attrs)?;
        }
    } else {
        for locus in loci {
            write_locus(out, locus)?;
            for transcript in annotation
                .transcripts
                .iter()
                .filter(|transcript| transcript.locus.as_deref() == Some(locus.id.as_str()))
            {
                write_transcript(out, transcript, attrs_filter, keep_all_attrs)?;
            }
        }
    }

    Ok(())
}

fn write_transcript<W: Write>(
    out: &mut W,
    transcript: &Transcript,
    attrs_filter: Option<&[String]>,
    keep_all_attrs: bool,
) -> io::Result<()> {
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

    if let Some(locus) = &transcript.locus {
        write!(out, ";locus={locus}")?;
    }

    if keep_all_attrs {
        for (attr_name, value) in &transcript.attrs {
            if should_emit_extra_attr(attr_name) && !value.is_empty() {
                write!(out, ";{attr_name}={value}")?;
            }
        }
    } else if let Some(attrs_filter) = attrs_filter {
        for attr_name in attrs_filter {
            if let Some(value) = transcript.attrs.get(attr_name) {
                if !value.is_empty() {
                    write!(out, ";{attr_name}={value}")?;
                }
            }
        }
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

fn should_emit_extra_attr(name: &str) -> bool {
    !matches!(
        name,
        "ID" | "Parent" | "geneID" | "gene_id" | "gene_name" | "transcript_id"
    )
}
