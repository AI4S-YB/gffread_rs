use std::collections::BTreeMap;
use std::fs;
use std::path::Path;

use crate::compat::CompatError;
use crate::model::{Annotation, Attrs, Gene, Segment, Transcript};
use crate::options::InputFormat;

pub fn load_annotation(path: &Path, input_format: InputFormat) -> Result<Annotation, CompatError> {
    let text = fs::read_to_string(path).map_err(|_| {
        CompatError::new(
            format!("Error: cannot open input file {}!\n", path.display()),
            1,
        )
    })?;
    if should_parse_bed(path, input_format) {
        parse_bed_annotation(&text)
    } else if should_parse_tlf(path, input_format) {
        parse_tlf_annotation(&text)
    } else {
        parse_annotation(&text)
    }
}

fn should_parse_tlf(path: &Path, input_format: InputFormat) -> bool {
    input_format == InputFormat::Tlf
        || path
            .extension()
            .and_then(|extension| extension.to_str())
            .is_some_and(|extension| extension.eq_ignore_ascii_case("tlf"))
}

fn should_parse_bed(path: &Path, input_format: InputFormat) -> bool {
    input_format == InputFormat::Bed
        || path
            .extension()
            .and_then(|extension| extension.to_str())
            .is_some_and(|extension| extension.eq_ignore_ascii_case("bed"))
}

pub fn parse_annotation(text: &str) -> Result<Annotation, CompatError> {
    let mut transcripts = Vec::<Transcript>::new();
    let mut genes = Vec::new();
    let mut gene_index = BTreeMap::<String, usize>::new();
    let mut transcript_index = BTreeMap::<String, usize>::new();
    let mut transcript_indices = BTreeMap::<String, Vec<usize>>::new();
    let mut ref_order = Vec::<String>::new();
    let mut header_comments = Vec::<String>::new();
    let mut seen_feature = false;
    let mut feature_order = 0usize;

    for line in text.lines() {
        if line.is_empty() {
            continue;
        }
        if line.starts_with('#') {
            if !seen_feature {
                header_comments.push(line.to_owned());
            }
            continue;
        }
        seen_feature = true;

        let fields: Vec<_> = line.split('\t').collect();
        if fields.len() != 9 {
            continue;
        }

        let attrs = parse_attrs(fields[8]);
        match fields[2] {
            "gene" => {
                if !ref_order.iter().any(|seqid| seqid == fields[0]) {
                    ref_order.push(fields[0].to_owned());
                }
                let id = attr(&attrs, "ID")
                    .or_else(|| attr(&attrs, "gene_id"))
                    .ok_or_else(|| CompatError::new("Error: gene record without ID\n", 1))?;
                let order = feature_order;
                feature_order += 1;
                genes.push(Gene {
                    seqid: fields[0].to_owned(),
                    source: fields[1].to_owned(),
                    feature: fields[2].to_owned(),
                    start: parse_u64(fields[3], "gene start")?,
                    end: parse_u64(fields[4], "gene end")?,
                    strand: fields[6].chars().next().unwrap_or('.'),
                    id,
                    order,
                    attrs,
                });
                let gene_idx = genes.len() - 1;
                if let Some(existing_idx) = gene_index.get(&genes[gene_idx].id).copied() {
                    if features_overlap(
                        genes[existing_idx].start,
                        genes[existing_idx].end,
                        genes[gene_idx].start,
                        genes[gene_idx].end,
                    ) {
                        eprintln!(
                            "Error: discarding overlapping duplicate gene feature ({}-{}) with ID={}",
                            genes[gene_idx].start, genes[gene_idx].end, genes[gene_idx].id
                        );
                    }
                } else {
                    gene_index.insert(genes[gene_idx].id.clone(), gene_idx);
                }
            }
            feature if is_transcript_feature(feature, &attrs) => {
                if !ref_order.iter().any(|seqid| seqid == fields[0]) {
                    ref_order.push(fields[0].to_owned());
                }

                let id = attr(&attrs, "ID")
                    .or_else(|| attr(&attrs, "transcript_id"))
                    .ok_or_else(|| CompatError::new("Error: transcript record without ID\n", 1))?;
                let order = feature_order;
                feature_order += 1;

                let transcript = Transcript {
                    seqid: fields[0].to_owned(),
                    source: fields[1].to_owned(),
                    feature: fields[2].to_owned(),
                    start: parse_u64(fields[3], "transcript start")?,
                    end: parse_u64(fields[4], "transcript end")?,
                    strand: fields[6].chars().next().unwrap_or('.'),
                    id: id.clone(),
                    order,
                    gene_id: attr(&attrs, "Parent")
                        .or_else(|| attr(&attrs, "geneID"))
                        .or_else(|| attr(&attrs, "gene_id")),
                    gene_name: attr(&attrs, "gene").or_else(|| attr(&attrs, "gene_name")),
                    attrs,
                    locus: None,
                    exons: Vec::new(),
                    cds: Vec::new(),
                };

                let transcript_idx = transcripts.len();
                if let Some(existing_idx) = transcript_index.get(&id).copied() {
                    if features_overlap(
                        transcripts[existing_idx].start,
                        transcripts[existing_idx].end,
                        transcript.start,
                        transcript.end,
                    ) {
                        eprintln!(
                            "Error: discarding overlapping duplicate {} feature ({}-{}) with ID={}",
                            feature, transcript.start, transcript.end, id
                        );
                    }
                } else {
                    transcript_index.insert(id.clone(), transcript_idx);
                }
                transcript_indices
                    .entry(id.clone())
                    .or_default()
                    .push(transcript_idx);
                transcripts.push(transcript);
            }
            "exon" | "CDS" => {
                let segment_start = parse_u64(fields[3], "segment start")?;
                let segment_end = parse_u64(fields[4], "segment end")?;
                let parent = attr(&attrs, "Parent")
                    .or_else(|| attr(&attrs, "transcript_id"))
                    .and_then(|value| value.split(',').next().map(str::to_owned))
                    .ok_or_else(|| CompatError::new("Error: child feature without Parent\n", 1))?;

                let Some(index) = find_parent_transcript(
                    &transcripts,
                    transcript_indices.get(&parent).map(Vec::as_slice),
                    fields[0],
                    fields[6].chars().next().unwrap_or('.'),
                    transcript_index.get(&parent).copied(),
                ) else {
                    continue;
                };

                let segment = Segment {
                    start: segment_start,
                    end: segment_end,
                    score: fields[5].to_owned(),
                    phase: fields[7].to_owned(),
                    attrs: filter_segment_attrs(&attrs),
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
        if transcript.gene_name.is_none() {
            if let Some(gene_id) = &transcript.gene_id {
                if let Some(gene) = gene_index.get(gene_id).and_then(|idx| genes.get(*idx)) {
                    transcript.gene_name = gene.gene_name().map(str::to_owned);
                }
            }
        }
        if transcript.exons.is_empty() {
            transcript.exons.push(Segment {
                start: transcript.start,
                end: transcript.end,
                score: ".".to_owned(),
                phase: ".".to_owned(),
                attrs: Attrs::new(),
            });
        }
        merge_adjacent_segments(&mut transcript.exons, transcript.strand, false);
        extend_exons_to_cds(&mut transcript.exons, &transcript.cds, transcript.strand);
        transcript
            .cds
            .sort_by_key(|segment| (segment.start, segment.end));
        dedup_exact_segments(&mut transcript.cds);
        if let Some(first_exon) = transcript.exons.first() {
            transcript.start = first_exon.start;
        }
        if let Some(last_exon) = transcript.exons.last() {
            transcript.end = last_exon.end;
        }
        update_missing_cds_phases(transcript);
    }

    for transcript in &mut transcripts {
        update_missing_cds_phases(transcript);
    }

    Ok(Annotation {
        transcripts,
        genes,
        ref_order,
        header_comments,
    })
}

fn features_overlap(left_start: u64, left_end: u64, right_start: u64, right_end: u64) -> bool {
    left_start <= right_end && right_start <= left_end
}

fn find_parent_transcript(
    transcripts: &[Transcript],
    candidates: Option<&[usize]>,
    child_seqid: &str,
    child_strand: char,
    fallback: Option<usize>,
) -> Option<usize> {
    let Some(candidates) = candidates else {
        return fallback;
    };

    candidates
        .iter()
        .copied()
        .find(|&index| {
            let transcript = &transcripts[index];
            transcript.seqid == child_seqid && strand_compatible(transcript.strand, child_strand)
        })
        .or(fallback)
}

fn strand_compatible(transcript_strand: char, child_strand: char) -> bool {
    transcript_strand == '.'
        || child_strand == '.'
        || transcript_strand == child_strand
        || !matches!(child_strand, '+' | '-')
}

fn is_transcript_feature(feature: &str, attrs: &Attrs) -> bool {
    matches!(
        feature,
        "mRNA"
            | "transcript"
            | "miRNA"
            | "ncRNA"
            | "tRNA"
            | "rRNA"
            | "snoRNA"
            | "snRNA"
            | "lnc_RNA"
            | "pseudogenic_transcript"
            | "primary_transcript"
    ) || (feature.ends_with("RNA") && attrs.get("ID").is_some() && has_parent_attr(attrs))
}

fn has_parent_attr(attrs: &Attrs) -> bool {
    attr(attrs, "Parent")
        .or_else(|| attr(attrs, "geneID"))
        .or_else(|| attr(attrs, "gene_id"))
        .is_some()
}

pub fn parse_bed_annotation(text: &str) -> Result<Annotation, CompatError> {
    let mut transcripts = Vec::new();
    let mut ref_order = Vec::<String>::new();
    let mut feature_order = 0usize;

    for line in text.lines() {
        if line.is_empty()
            || line.starts_with('#')
            || line.starts_with("browser ")
            || line.starts_with("track ")
        {
            continue;
        }

        let fields: Vec<_> = line.split('\t').collect();
        if fields.len() < 3 {
            continue;
        }

        let Some(mut transcript) = parse_bed_transcript(&fields)? else {
            continue;
        };
        transcript.order = feature_order;
        feature_order += 1;

        if !ref_order.iter().any(|seqid| seqid == &transcript.seqid) {
            ref_order.push(transcript.seqid.clone());
        }
        transcripts.push(transcript);
    }

    Ok(Annotation {
        transcripts,
        genes: Vec::new(),
        ref_order,
        header_comments: Vec::new(),
    })
}

pub fn parse_tlf_annotation(text: &str) -> Result<Annotation, CompatError> {
    let mut annotation = parse_annotation(text)?;
    for transcript in &mut annotation.transcripts {
        if let Some(exons) = transcript.attrs.get("exons") {
            transcript.exons = parse_segment_list(exons)?;
        }
        if let Some(cds) = transcript.attrs.get("CDS") {
            transcript.cds = parse_tlf_cds_segments(
                cds,
                &transcript.exons,
                transcript
                    .attrs
                    .get("CDSphase")
                    .map(str::to_owned)
                    .unwrap_or_else(|| "0".to_owned()),
            )?;
        }
        transcript
            .cds
            .sort_by_key(|segment| (segment.start, segment.end));
        dedup_exact_segments(&mut transcript.cds);
        merge_adjacent_segments(&mut transcript.exons, transcript.strand, false);
        if let Some(first_exon) = transcript.exons.first() {
            transcript.start = first_exon.start;
        }
        if let Some(last_exon) = transcript.exons.last() {
            transcript.end = last_exon.end;
        }
        update_missing_cds_phases(transcript);
    }

    Ok(annotation)
}

fn merge_adjacent_segments(segments: &mut Vec<Segment>, strand: char, keep_phase: bool) {
    segments.sort_by_key(|segment| (segment.start, segment.end));
    let mut merged: Vec<Segment> = Vec::with_capacity(segments.len());
    for mut segment in segments.drain(..) {
        if !keep_phase {
            segment.phase = ".".to_owned();
        }
        if let Some(last) = merged.last_mut() {
            if segment.start <= last.end + 1 {
                if keep_phase {
                    if segment.start < last.start && strand == '+' {
                        last.phase = segment.phase.clone();
                    }
                    if segment.end > last.end && strand == '-' {
                        last.phase = segment.phase.clone();
                    }
                }
                if segment.start < last.start {
                    last.start = segment.start;
                }
                if segment.end > last.end {
                    last.end = segment.end;
                }
                continue;
            }
        }
        merged.push(segment);
    }
    *segments = merged;
}

fn dedup_exact_segments(segments: &mut Vec<Segment>) {
    let mut deduped = Vec::with_capacity(segments.len());
    for segment in segments.drain(..) {
        if deduped
            .last()
            .is_some_and(|last: &Segment| last.start == segment.start && last.end == segment.end)
        {
            continue;
        }
        deduped.push(segment);
    }
    *segments = deduped;
}

fn extend_exons_to_cds(exons: &mut Vec<Segment>, cds_segments: &[Segment], strand: char) {
    for cds in cds_segments {
        if let Some(exon) = exons
            .iter_mut()
            .find(|exon| features_overlap(exon.start, exon.end, cds.start, cds.end))
        {
            if cds.start < exon.start {
                exon.start = cds.start;
            }
            if cds.end > exon.end {
                exon.end = cds.end;
            }
        } else {
            exons.push(Segment {
                start: cds.start,
                end: cds.end,
                score: cds.score.clone(),
                phase: ".".to_owned(),
                attrs: Attrs::new(),
            });
        }
    }
    merge_adjacent_segments(exons, strand, false);
}

fn update_missing_cds_phases(transcript: &mut Transcript) {
    if transcript.cds.is_empty() {
        return;
    }

    let initial_phase = if transcript.strand == '-' {
        transcript.cds.last()
    } else {
        transcript.cds.first()
    }
    .and_then(|segment| segment.phase.chars().next());

    let has_invalid_phase = transcript.cds.iter().any(|segment| {
        !segment
            .phase
            .chars()
            .next()
            .is_some_and(|phase| matches!(phase, '0' | '1' | '2'))
    });
    if !has_invalid_phase && !cds_exon_compatible(transcript) {
        return;
    }

    let cdsacc = match initial_phase {
        Some('1') => 2,
        Some('2') => 1,
        _ => 0,
    };
    recompute_cds_phases(transcript, cdsacc);
}

fn cds_exon_compatible(transcript: &Transcript) -> bool {
    if transcript.exons.is_empty() || transcript.cds.is_empty() {
        return false;
    }
    if transcript.cds.len() == 1 {
        let cds = &transcript.cds[0];
        return transcript
            .exons
            .iter()
            .any(|exon| cds.start >= exon.start && cds.end <= exon.end);
    }

    let Some(mut exon_index) = transcript.exons.iter().position(|exon| {
        transcript.cds[0].start >= exon.start && transcript.cds[0].start <= exon.end
    }) else {
        return false;
    };

    let mut cds_index = 0usize;
    while exon_index + 1 < transcript.exons.len() && cds_index + 1 < transcript.cds.len() {
        if transcript.exons[exon_index].end != transcript.cds[cds_index].end
            || transcript.exons[exon_index + 1].start != transcript.cds[cds_index + 1].start
        {
            return false;
        }
        exon_index += 1;
        cds_index += 1;
    }

    let last_cds = transcript.cds.last().expect("checked non-empty");
    cds_index + 1 == transcript.cds.len()
        && last_cds.end >= transcript.exons[exon_index].start
        && last_cds.end <= transcript.exons[exon_index].end
}

fn recompute_cds_phases(transcript: &mut Transcript, mut cdsacc: u64) {
    if transcript.strand == '-' {
        for index in (0..transcript.cds.len()).rev() {
            let phase = (3 - (cdsacc % 3)) % 3;
            transcript.cds[index].phase = phase.to_string();
            cdsacc += transcript.cds[index].end - transcript.cds[index].start + 1;
        }
    } else {
        for segment in &mut transcript.cds {
            let phase = (3 - (cdsacc % 3)) % 3;
            segment.phase = phase.to_string();
            cdsacc += segment.end - segment.start + 1;
        }
    }
}

fn parse_segment_list(raw: &str) -> Result<Vec<Segment>, CompatError> {
    let mut segments = Vec::new();
    for part in raw.split(',') {
        let part = part.trim();
        if part.is_empty() {
            continue;
        }
        let Some((start, end)) = part.split_once('-') else {
            return Ok(Vec::new());
        };
        segments.push(Segment {
            start: parse_u64(start, "segment start")?,
            end: parse_u64(end, "segment end")?,
            score: ".".to_owned(),
            phase: ".".to_owned(),
            attrs: Attrs::new(),
        });
    }
    segments.sort_by_key(|segment| (segment.start, segment.end));
    Ok(segments)
}

fn parse_tlf_cds_segments(
    raw: &str,
    exons: &[Segment],
    phase: String,
) -> Result<Vec<Segment>, CompatError> {
    if let Some((start, end)) = raw.split_once(':') {
        let cds_start = parse_u64(start, "TLF CDS start")?;
        let cds_end = parse_u64(end, "TLF CDS end")?;
        return Ok(bed_cds_segments(exons, cds_start, cds_end, phase));
    }

    let mut segments = parse_segment_list(raw)?;
    if let Some(first) = segments.first_mut() {
        first.phase = phase;
    }
    Ok(segments)
}

fn parse_bed_transcript(fields: &[&str]) -> Result<Option<Transcript>, CompatError> {
    let start0 = parse_u64(fields[1], "BED start")?;
    let end = parse_u64(fields[2], "BED end")?;
    let (start, end) = if end < start0 + 1 {
        (end, start0 + 1)
    } else {
        (start0 + 1, end)
    };
    let id = fields
        .get(3)
        .filter(|value| !value.is_empty())
        .unwrap_or(&fields[0]);
    let strand = fields
        .get(5)
        .and_then(|value| value.chars().next())
        .filter(|value| matches!(value, '+' | '-' | '.'))
        .unwrap_or('.');
    let attrs = fields
        .get(12)
        .map_or_else(Attrs::new, |raw| parse_attrs(raw));

    let exons = if fields.len() > 11 {
        let block_count = parse_u64(fields[9], "BED block count")? as usize;
        let sizes = parse_bed_blocks(fields[10]);
        let starts = parse_bed_blocks(fields[11]);
        if block_count == 0 || sizes.len() < block_count || starts.len() < block_count {
            return Ok(None);
        }
        let mut exons = Vec::with_capacity(block_count);
        for index in 0..block_count {
            let exon_start = start + starts[index];
            let exon_end = exon_start + sizes[index].saturating_sub(1);
            if sizes[index] == 0 || exon_start > end || exon_end > end {
                return Ok(None);
            }
            exons.push(Segment {
                start: exon_start,
                end: exon_end,
                score: ".".to_owned(),
                phase: ".".to_owned(),
                attrs: Attrs::new(),
            });
        }
        exons
    } else {
        vec![Segment {
            start,
            end,
            score: ".".to_owned(),
            phase: ".".to_owned(),
            attrs: Attrs::new(),
        }]
    };

    let (cds_start, cds_end, cds_phase) = bed_cds_span(fields, &attrs, start, end)?;
    let cds = match (cds_start, cds_end) {
        (Some(cds_start), Some(cds_end)) => bed_cds_segments(&exons, cds_start, cds_end, cds_phase),
        _ => Vec::new(),
    };

    Ok(Some(Transcript {
        seqid: fields[0].to_owned(),
        source: "BED".to_owned(),
        feature: "transcript".to_owned(),
        start,
        end,
        strand,
        id: (*id).to_owned(),
        order: 0,
        gene_id: None,
        gene_name: None,
        attrs,
        locus: None,
        exons,
        cds,
    }))
}

fn parse_bed_blocks(raw: &str) -> Vec<u64> {
    raw.split(',')
        .filter(|part| !part.is_empty())
        .filter_map(|part| part.parse::<u64>().ok())
        .collect()
}

fn bed_cds_span(
    fields: &[&str],
    attrs: &Attrs,
    start: u64,
    end: u64,
) -> Result<(Option<u64>, Option<u64>, String), CompatError> {
    let mut phase = attr(attrs, "CDSphase").unwrap_or_else(|| "0".to_owned());
    if !matches!(phase.as_str(), "0" | "1" | "2") {
        phase = "0".to_owned();
    }

    if let Some(cds_attr) = attr(attrs, "CDS") {
        if let Some((raw_start, raw_end)) = cds_attr.split_once(':') {
            let cds_start0 = parse_u64(raw_start, "BED CDS start")?;
            let cds_end = parse_u64(raw_end, "BED CDS end")?;
            let cds_start = cds_start0 + 1;
            if cds_start >= start && cds_end <= end && cds_end >= cds_start {
                return Ok((Some(cds_start), Some(cds_end), phase));
            }
        }
    }

    if fields.len() > 7 {
        let cds_start0 = parse_u64(fields[6], "BED thick start")?;
        let cds_end = parse_u64(fields[7], "BED thick end")?;
        let cds_start = cds_start0 + 1;
        if cds_end > cds_start0 && cds_start >= start && cds_end <= end {
            return Ok((Some(cds_start), Some(cds_end), phase));
        }
    }

    Ok((None, None, phase))
}

fn bed_cds_segments(
    exons: &[Segment],
    cds_start: u64,
    cds_end: u64,
    first_phase: String,
) -> Vec<Segment> {
    let mut cds_segments = Vec::new();
    let mut offset = 0u64;
    for exon in exons {
        let start = exon.start.max(cds_start);
        let end = exon.end.min(cds_end);
        if start <= end {
            let phase = if cds_segments.is_empty() {
                first_phase.clone()
            } else {
                ((3 - (offset % 3)) % 3).to_string()
            };
            offset += end - start + 1;
            cds_segments.push(Segment {
                start,
                end,
                score: ".".to_owned(),
                phase,
                attrs: Attrs::new(),
            });
        }
    }
    cds_segments
}

fn parse_u64(value: &str, field_name: &str) -> Result<u64, CompatError> {
    value
        .parse()
        .map_err(|_| CompatError::new(format!("Error: invalid {field_name}\n"), 1))
}

fn parse_attrs(raw: &str) -> Attrs {
    let mut attrs = Attrs::new();

    for part in raw.split(';') {
        let part = part.trim();
        if part.is_empty() {
            continue;
        }

        if let Some((key, value)) = part.split_once('=') {
            attrs.insert_or_replace(
                key.trim().to_owned(),
                value.trim().trim_matches('"').to_owned(),
            );
            continue;
        }

        if let Some((key, value)) = part.split_once(' ') {
            attrs.insert_or_replace(
                key.trim().to_owned(),
                value.trim().trim_matches('"').to_owned(),
            );
        }
    }

    attrs
}

fn attr(attrs: &Attrs, name: &str) -> Option<String> {
    attrs.get(name).map(str::to_owned)
}

fn filter_segment_attrs(attrs: &Attrs) -> Attrs {
    let mut filtered = Attrs::new();
    for attr in attrs.iter() {
        if matches!(attr.key.as_str(), "Parent" | "transcript_id") {
            continue;
        }
        filtered.insert_or_replace(attr.key.clone(), attr.value.clone());
    }
    filtered
}
