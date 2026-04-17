use std::collections::BTreeMap;
use std::fs;
use std::path::Path;

use crate::compat::CompatError;
use crate::model::{Annotation, Segment, Transcript};

pub fn load_annotation(path: &Path) -> Result<Annotation, CompatError> {
    let text = fs::read_to_string(path).map_err(|_| {
        CompatError::new(
            format!("Error: cannot open input file {}!\n", path.display()),
            1,
        )
    })?;
    parse_annotation(&text)
}

pub fn parse_annotation(text: &str) -> Result<Annotation, CompatError> {
    let mut transcripts = Vec::new();
    let mut transcript_index = BTreeMap::<String, usize>::new();

    for line in text.lines() {
        if line.is_empty() || line.starts_with('#') {
            continue;
        }

        let fields: Vec<_> = line.split('\t').collect();
        if fields.len() != 9 {
            continue;
        }

        let attrs = parse_attrs(fields[8]);
        match fields[2] {
            "mRNA" | "transcript" => {
                let id = attr(&attrs, "ID")
                    .or_else(|| attr(&attrs, "transcript_id"))
                    .ok_or_else(|| CompatError::new("Error: transcript record without ID\n", 1))?;

                let transcript = Transcript {
                    seqid: fields[0].to_owned(),
                    source: fields[1].to_owned(),
                    feature: fields[2].to_owned(),
                    start: parse_u64(fields[3], "transcript start")?,
                    end: parse_u64(fields[4], "transcript end")?,
                    strand: fields[6].chars().next().unwrap_or('.'),
                    id: id.clone(),
                    gene_id: attr(&attrs, "Parent")
                        .or_else(|| attr(&attrs, "geneID"))
                        .or_else(|| attr(&attrs, "gene_id")),
                    gene_name: attr(&attrs, "gene")
                        .or_else(|| attr(&attrs, "gene_name"))
                        .or_else(|| attr(&attrs, "Name")),
                    attrs,
                    exons: Vec::new(),
                    cds: Vec::new(),
                };

                transcript_index.insert(id, transcripts.len());
                transcripts.push(transcript);
            }
            "exon" | "CDS" => {
                let parent = attr(&attrs, "Parent")
                    .or_else(|| attr(&attrs, "transcript_id"))
                    .and_then(|value| value.split(',').next().map(str::to_owned))
                    .ok_or_else(|| CompatError::new("Error: child feature without Parent\n", 1))?;

                let Some(index) = transcript_index.get(&parent).copied() else {
                    continue;
                };

                let segment = Segment {
                    start: parse_u64(fields[3], "segment start")?,
                    end: parse_u64(fields[4], "segment end")?,
                    phase: fields[7].to_owned(),
                };

                if fields[2] == "exon" {
                    transcripts[index].exons.push(segment);
                } else {
                    transcripts[index].cds.push(segment);
                }
            }
            _ => {}
        }
    }

    for transcript in &mut transcripts {
        transcript
            .exons
            .sort_by_key(|segment| (segment.start, segment.end));
        transcript
            .cds
            .sort_by_key(|segment| (segment.start, segment.end));
    }

    transcripts.sort_by(|left, right| {
        (&left.seqid, left.start, left.end, &left.id).cmp(&(
            &right.seqid,
            right.start,
            right.end,
            &right.id,
        ))
    });

    Ok(Annotation { transcripts })
}

fn parse_u64(value: &str, field_name: &str) -> Result<u64, CompatError> {
    value
        .parse()
        .map_err(|_| CompatError::new(format!("Error: invalid {field_name}\n"), 1))
}

fn parse_attrs(raw: &str) -> BTreeMap<String, String> {
    let mut attrs = BTreeMap::new();

    for part in raw.split(';') {
        let part = part.trim();
        if part.is_empty() {
            continue;
        }

        if let Some((key, value)) = part.split_once('=') {
            attrs.insert(
                key.trim().to_owned(),
                value.trim().trim_matches('"').to_owned(),
            );
            continue;
        }

        if let Some((key, value)) = part.split_once(' ') {
            attrs.insert(
                key.trim().to_owned(),
                value.trim().trim_matches('"').to_owned(),
            );
        }
    }

    attrs
}

fn attr(attrs: &BTreeMap<String, String>, name: &str) -> Option<String> {
    attrs.get(name).filter(|value| !value.is_empty()).cloned()
}
