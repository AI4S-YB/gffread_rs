use std::io::{self, Write};

use crate::emit::locus::write_locus;
use crate::model::{Annotation, Attrs, Gene, Locus, Segment, Transcript};

pub struct Gff3Options<'a> {
    pub version: &'a str,
    pub command_line: &'a str,
    pub track_label: Option<&'a str>,
    pub attrs_filter: Option<&'a [String]>,
    pub keep_all_attrs: bool,
    pub gather_exon_attrs: bool,
    pub keep_exon_attrs: bool,
    pub keep_genes: bool,
    pub keep_comments: bool,
    pub decode_attrs: bool,
}

pub fn write_gff3<W: Write>(
    out: &mut W,
    annotation: &Annotation,
    loci: &[Locus],
    options: &Gff3Options<'_>,
) -> io::Result<()> {
    if options.keep_comments {
        for comment in &annotation.header_comments {
            writeln!(out, "{comment}")?;
        }
    } else {
        writeln!(out, "##gff-version 3")?;
        writeln!(out, "# gffread v{}", options.version)?;
        writeln!(out, "# {}", options.command_line)?;
    }

    if loci.is_empty() {
        write_records(out, annotation, None, options)?;
    } else {
        for locus in loci {
            write_locus(out, locus)?;
            write_records(out, annotation, Some(locus.id.as_str()), options)?;
        }
    }

    Ok(())
}

fn write_records<W: Write>(
    out: &mut W,
    annotation: &Annotation,
    locus_id: Option<&str>,
    options: &Gff3Options<'_>,
) -> io::Result<()> {
    let transcripts = annotation
        .transcripts
        .iter()
        .filter(|transcript| transcript.locus.as_deref() == locus_id);

    let mut emitted_genes = Vec::<String>::new();
    for transcript in transcripts {
        if options.keep_genes {
            if let Some(gene_id) = transcript.gene_id.as_deref() {
                if !emitted_genes.iter().any(|id| id == gene_id) {
                    if let Some(gene) = annotation.genes.iter().find(|gene| gene.id == gene_id) {
                        write_gene(
                            out,
                            gene,
                            options.track_label,
                            options.keep_all_attrs,
                            options.decode_attrs,
                        )?;
                        emitted_genes.push(gene_id.to_owned());
                    }
                }
            }
        }
        write_transcript(out, transcript, options)?;
    }

    Ok(())
}

fn write_transcript<W: Write>(
    out: &mut W,
    transcript: &Transcript,
    options: &Gff3Options<'_>,
) -> io::Result<()> {
    let source = options.track_label.unwrap_or(&transcript.source);
    let transcript_attrs = merged_transcript_attrs(
        transcript,
        options.keep_all_attrs,
        options.gather_exon_attrs,
        options.keep_exon_attrs,
    );
    write!(
        out,
        "{}\t{}\t{}\t{}\t{}\t.\t{}\t.\tID={}",
        transcript.seqid,
        source,
        transcript.feature,
        transcript.start,
        transcript.end,
        transcript.strand,
        transcript.id
    )?;

    if let Some(gene_id) = &transcript.gene_id {
        if options.keep_genes {
            write!(out, ";Parent={}", attr_value(gene_id, options.decode_attrs))?;
        } else {
            write!(out, ";geneID={}", attr_value(gene_id, options.decode_attrs))?;
        }
    }

    if let Some(gene_name) = &transcript.gene_name {
        write!(
            out,
            ";gene_name={}",
            attr_value(gene_name, options.decode_attrs)
        )?;
    }

    if let Some(locus) = &transcript.locus {
        write!(out, ";locus={locus}")?;
    }

    if options.keep_all_attrs {
        for attr in transcript_attrs.iter() {
            if should_emit_extra_attr(&attr.key) && !attr.value.is_empty() {
                write!(
                    out,
                    ";{}={}",
                    attr.key,
                    attr_value(&attr.value, options.decode_attrs)
                )?;
            }
        }
    } else if let Some(attrs_filter) = options.attrs_filter {
        for attr_name in attrs_filter {
            if let Some(value) = transcript_attrs.get(attr_name) {
                if !value.is_empty() {
                    write!(
                        out,
                        ";{attr_name}={}",
                        attr_value(value, options.decode_attrs)
                    )?;
                }
            }
        }
    }

    writeln!(out)?;

    for exon in &transcript.exons {
        write_segment(
            out,
            transcript,
            source,
            exon,
            &transcript.exons,
            "exon",
            options,
        )?;
    }

    for cds in &transcript.cds {
        write_segment(
            out,
            transcript,
            source,
            cds,
            &transcript.cds,
            "CDS",
            options,
        )?;
    }

    Ok(())
}

fn write_segment<W: Write>(
    out: &mut W,
    transcript: &Transcript,
    source: &str,
    segment: &Segment,
    sibling_segments: &[Segment],
    feature: &str,
    options: &Gff3Options<'_>,
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
    write!(
        out,
        "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\tParent={}",
        transcript.seqid,
        source,
        feature,
        segment.start,
        segment.end,
        score,
        transcript.strand,
        phase,
        transcript.id
    )?;
    if options.keep_all_attrs {
        for attr in segment.attrs.iter() {
            if !options.keep_exon_attrs
                && should_skip_segment_attr(
                    attr.key.as_str(),
                    attr.value.as_str(),
                    sibling_segments,
                    options.gather_exon_attrs,
                )
            {
                continue;
            }
            if !attr.value.is_empty() {
                write!(
                    out,
                    ";{}={}",
                    attr.key,
                    attr_value(&attr.value, options.decode_attrs)
                )?;
            }
        }
    }
    writeln!(out)
}

fn write_gene<W: Write>(
    out: &mut W,
    gene: &Gene,
    track_label: Option<&str>,
    keep_all_attrs: bool,
    decode_attrs: bool,
) -> io::Result<()> {
    let source = track_label.unwrap_or(&gene.source);
    write!(
        out,
        "{}\t{}\t{}\t{}\t{}\t.\t{}\t.\tID={}",
        gene.seqid, source, gene.feature, gene.start, gene.end, gene.strand, gene.id
    )?;
    if keep_all_attrs {
        for attr in gene.attrs.iter() {
            if should_emit_gene_attr(attr.key.as_str()) && !attr.value.is_empty() {
                write!(
                    out,
                    ";{}={}",
                    attr.key,
                    attr_value(&attr.value, decode_attrs)
                )?;
            }
        }
    }
    writeln!(out)
}

fn should_emit_extra_attr(name: &str) -> bool {
    !matches!(
        name,
        "ID" | "Parent" | "geneID" | "gene_id" | "gene_name" | "transcript_id"
    )
}

fn should_emit_gene_attr(name: &str) -> bool {
    !matches!(name, "ID" | "geneID" | "gene_id" | "gene_name")
}

fn merged_transcript_attrs(
    transcript: &Transcript,
    keep_all_attrs: bool,
    gather_exon_attrs: bool,
    keep_exon_attrs: bool,
) -> Attrs {
    let mut attrs = transcript.attrs.clone();
    if !keep_all_attrs || keep_exon_attrs {
        return attrs;
    }

    if gather_exon_attrs {
        copy_first_segment_attrs(&mut attrs, &transcript.exons);
        copy_first_segment_attrs(&mut attrs, &transcript.cds);
    } else {
        reduce_segment_attrs(&mut attrs, &transcript.exons, "exon_");
        reduce_segment_attrs(&mut attrs, &transcript.cds, "CDS_");
    }

    attrs
}

fn copy_first_segment_attrs(attrs: &mut Attrs, segments: &[Segment]) {
    let Some(first) = segments.first() else {
        return;
    };
    for attr in first.attrs.iter() {
        if attr.key.starts_with("exon_") || attr.key == "exon" {
            continue;
        }
        attrs.push_unique(attr.key.clone(), attr.value.clone());
    }
}

fn reduce_segment_attrs(attrs: &mut Attrs, segments: &[Segment], prefix: &str) {
    let Some(first) = segments.first() else {
        return;
    };
    for attr in first.attrs.iter() {
        let discard_all = matches!(attr.key.as_str(), "exon_id" | "exon_number");
        let same_in_all = !discard_all
            && segments
                .iter()
                .skip(1)
                .all(|segment| segment.attrs.get(&attr.key) == Some(attr.value.as_str()));
        if !same_in_all {
            continue;
        }

        if let Some(existing) = attrs.get(&attr.key) {
            if existing != attr.value {
                attrs.push_unique(format!("{prefix}{}", attr.key), attr.value.clone());
            }
        } else {
            attrs.push_unique(attr.key.clone(), attr.value.clone());
        }
    }
}

fn should_skip_segment_attr(
    key: &str,
    value: &str,
    sibling_segments: &[Segment],
    gather_exon_attrs: bool,
) -> bool {
    if gather_exon_attrs {
        return true;
    }
    if matches!(key, "exon_id" | "exon_number") {
        return true;
    }
    let Some(first) = sibling_segments.first() else {
        return false;
    };
    first
        .attrs
        .get(key)
        .is_some_and(|first_value| first_value == value)
        && sibling_segments
            .iter()
            .skip(1)
            .all(|segment| segment.attrs.get(key) == Some(value))
}

fn attr_value(value: &str, decode_attrs: bool) -> String {
    if decode_attrs {
        decode_hex_chars(value)
    } else {
        value.to_owned()
    }
}

fn decode_hex_chars(value: &str) -> String {
    let mut decoded = String::new();
    let bytes = value.as_bytes();
    let mut index = 0usize;
    while index < bytes.len() {
        if bytes[index] == b'%' && index + 2 < bytes.len() {
            if let (Some(high), Some(low)) =
                (hex_value(bytes[index + 1]), hex_value(bytes[index + 2]))
            {
                let mut ch = (high << 4) + low;
                index += 3;
                if ch == b'%' {
                    decoded.push_str("prc");
                    continue;
                }
                if ch == b';' {
                    ch = b'.';
                } else if ch <= b'\t' {
                    ch = b' ';
                }
                if ch >= b' ' {
                    decoded.push(ch as char);
                    continue;
                }
                decoded.push(bytes[index - 1] as char);
                continue;
            }
        }
        decoded.push(bytes[index] as char);
        index += 1;
    }
    decoded
}

fn hex_value(byte: u8) -> Option<u8> {
    match byte {
        b'0'..=b'9' => Some(byte - b'0'),
        b'a'..=b'f' => Some(byte - b'a' + 10),
        b'A'..=b'F' => Some(byte - b'A' + 10),
        _ => None,
    }
}
