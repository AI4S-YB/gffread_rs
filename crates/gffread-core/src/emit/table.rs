use std::io::{self, Write};

use crate::model::{Annotation, Transcript};

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum TableField {
    Id,
    GeneId,
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
        .map(|part| {
            if let Some(special) = part.strip_prefix('@') {
                // Special @fields are case-insensitive, like the C++ parser
                // which lowercases the token before the lookup.
                match special.to_ascii_lowercase().as_str() {
                    "id" => return TableField::Id,
                    "geneid" => return TableField::GeneId,
                    "chr" => return TableField::Chr,
                    "start" => return TableField::Start,
                    "end" => return TableField::End,
                    "strand" => return TableField::Strand,
                    "exons" => return TableField::Exons,
                    _ => {}
                }
            }
            match part {
                "ID" | "transcript_id" => TableField::Id,
                "geneID" | "gene_id" => TableField::GeneId,
                other => TableField::Attr(other.to_owned()),
            }
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
        TableField::GeneId => transcript.gene_id.as_deref().unwrap_or(".").to_owned(),
        TableField::Chr => transcript.seqid.clone(),
        TableField::Start => transcript.start.to_string(),
        TableField::End => transcript.end.to_string(),
        TableField::Strand => transcript.strand.to_string(),
        TableField::Exons => transcript.exon_list(),
        TableField::Attr(name) => transcript.attrs.get(name).unwrap_or_default().to_owned(),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::model::{Annotation, Attrs, Segment};

    fn transcript(id: &str, gene_id: Option<&str>) -> Transcript {
        Transcript {
            seqid: "chr1".to_owned(),
            source: "test".to_owned(),
            feature: "mRNA".to_owned(),
            start: 1,
            end: 100,
            strand: '+',
            id: id.to_owned(),
            order: 0,
            gene_id: gene_id.map(str::to_owned),
            gene_name: None,
            attrs: Attrs::new(),
            locus: None,
            exons: vec![Segment {
                start: 1,
                end: 100,
                score: ".".to_owned(),
                phase: ".".to_owned(),
                attrs: Attrs::new(),
            }],
            cds: Vec::new(),
        }
    }

    #[test]
    fn parses_gene_id_special_field_case_insensitively() {
        assert_eq!(parse_format("@geneid"), vec![TableField::GeneId]);
        assert_eq!(parse_format("@GENEID"), vec![TableField::GeneId]);
        assert_eq!(parse_format("@GeneId"), vec![TableField::GeneId]);
    }

    #[test]
    fn parses_gene_id_attribute_aliases() {
        assert_eq!(parse_format("geneID"), vec![TableField::GeneId]);
        assert_eq!(parse_format("gene_id"), vec![TableField::GeneId]);
    }

    #[test]
    fn writes_gene_id_and_dot_placeholder_when_missing() {
        let annotation = Annotation {
            transcripts: vec![transcript("rna1", Some("gene1")), transcript("rna2", None)],
            genes: Vec::new(),
            ref_order: Vec::new(),
            header_comments: Vec::new(),
        };
        let mut out = Vec::new();
        write_table(&mut out, &annotation, "@id,@GENEID").unwrap();
        assert_eq!(String::from_utf8(out).unwrap(), "rna1\tgene1\nrna2\t.\n");
    }
}
