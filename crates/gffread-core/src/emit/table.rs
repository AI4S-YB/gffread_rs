use std::io::{self, Write};

use crate::model::{Annotation, Transcript};

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum TableField {
    Id,
    Chr,
    Start,
    End,
    Strand,
    Exons,
    Attr(String),
}

pub fn parse_format(format: &str) -> Vec<TableField> {
    format
        .split([',', ';', ':', ' '])
        .filter(|part| !part.is_empty())
        .map(|part| match part {
            "@id" | "ID" | "transcript_id" => TableField::Id,
            "@chr" => TableField::Chr,
            "@start" => TableField::Start,
            "@end" => TableField::End,
            "@strand" => TableField::Strand,
            "@exons" => TableField::Exons,
            other => TableField::Attr(other.to_owned()),
        })
        .collect()
}

pub fn write_table<W: Write>(out: &mut W, annotation: &Annotation, format: &str) -> io::Result<()> {
    let fields = parse_format(format);

    for transcript in &annotation.transcripts {
        let row = fields
            .iter()
            .map(|field| value_for(transcript, field))
            .collect::<Vec<_>>()
            .join("\t");
        writeln!(out, "{row}")?;
    }

    Ok(())
}

fn value_for(transcript: &Transcript, field: &TableField) -> String {
    match field {
        TableField::Id => transcript.id.clone(),
        TableField::Chr => transcript.seqid.clone(),
        TableField::Start => transcript.start.to_string(),
        TableField::End => transcript.end.to_string(),
        TableField::Strand => transcript.strand.to_string(),
        TableField::Exons => transcript.exon_list(),
        TableField::Attr(name) => transcript.attrs.get(name).unwrap_or_default().to_owned(),
    }
}
