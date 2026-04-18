use std::collections::BTreeMap;
use std::fmt::Write as _;
use std::fs;
use std::io::{self, Write};
use std::path::{Path, PathBuf};

use crate::compat::CompatError;
use crate::model::{Annotation, Segment, Transcript};

pub type Genome = BTreeMap<String, Vec<u8>>;

pub fn load_genome(path: &Path) -> Result<Genome, CompatError> {
    let bytes = fs::read(path).map_err(|_| {
        CompatError::new(
            format!(
                "Error: couldn't open genomic sequences file {}\n",
                path.display()
            ),
            1,
        )
    })?;

    let mut genome = Genome::new();
    let mut index_entries = Vec::new();
    let mut current_id = None::<String>;
    let mut current_seq = Vec::new();
    let mut current_offset = 0usize;
    let mut current_line_bases = 0usize;
    let mut current_line_width = 0usize;
    let mut position = 0usize;

    for raw_line in bytes.split_inclusive(|byte| *byte == b'\n') {
        let line = trim_line_ending(raw_line);
        if let Some(rest) = line.strip_prefix(b">") {
            finalize_record(
                &mut genome,
                &mut index_entries,
                &mut current_id,
                &mut current_seq,
                current_offset,
                current_line_bases,
                current_line_width,
            );

            current_id = Some(parse_record_id(rest));
            current_offset = position + raw_line.len();
            current_line_bases = 0;
            current_line_width = 0;
        } else if current_id.is_some() {
            let bases = line
                .iter()
                .copied()
                .filter(|byte| !byte.is_ascii_whitespace())
                .map(|byte| byte.to_ascii_uppercase())
                .collect::<Vec<_>>();

            if !bases.is_empty() {
                if current_line_bases == 0 {
                    current_line_bases = bases.len();
                    current_line_width = raw_line.len();
                }
                current_seq.extend_from_slice(&bases);
            }
        }

        position += raw_line.len();
    }

    finalize_record(
        &mut genome,
        &mut index_entries,
        &mut current_id,
        &mut current_seq,
        current_offset,
        current_line_bases,
        current_line_width,
    );

    write_fasta_index(path, &index_entries)?;
    Ok(genome)
}

pub fn write_transcript_fasta<W: Write>(
    out: &mut W,
    annotation: &Annotation,
    genome: &Genome,
    padding: u64,
    suppress_cds: bool,
    write_segments: bool,
) -> Result<(), CompatError> {
    for transcript in &annotation.transcripts {
        let seq = spliced_sequence_with_padding(transcript, &transcript.exons, genome, padding)?;
        let defline =
            transcript_defline(transcript, genome, padding, suppress_cds, write_segments)?;
        write_fasta_record(out, &defline, &seq, false)
            .map_err(|err| CompatError::new(format!("Error writing FASTA: {err}\n"), 1))?;
    }
    Ok(())
}

pub fn write_unspliced_fasta<W: Write>(
    out: &mut W,
    annotation: &Annotation,
    genome: &Genome,
    padding: u64,
) -> Result<(), CompatError> {
    for transcript in &annotation.transcripts {
        let seq = unspliced_sequence(transcript, genome, padding)?;
        write_fasta_record(out, &transcript.id, &seq, false)
            .map_err(|err| CompatError::new(format!("Error writing FASTA: {err}\n"), 1))?;
    }
    Ok(())
}

pub fn write_cds_fasta<W: Write>(
    out: &mut W,
    annotation: &Annotation,
    genome: &Genome,
    write_segments: bool,
) -> Result<(), CompatError> {
    for transcript in &annotation.transcripts {
        if transcript.cds.is_empty() {
            continue;
        }

        let seq = spliced_cds_sequence(transcript, genome)?;
        let defline = if write_segments {
            projected_defline(transcript, &transcript.cds)
        } else {
            transcript.id.clone()
        };
        write_fasta_record(out, &defline, &seq, false)
            .map_err(|err| CompatError::new(format!("Error writing FASTA: {err}\n"), 1))?;
    }
    Ok(())
}

pub fn write_protein_fasta<W: Write>(
    out: &mut W,
    annotation: &Annotation,
    genome: &Genome,
    use_star_stop: bool,
    write_segments: bool,
) -> Result<(), CompatError> {
    for transcript in &annotation.transcripts {
        if transcript.cds.is_empty() {
            continue;
        }

        let cds = spliced_cds_sequence(transcript, genome)?;
        let protein = translate(&cds);
        if protein.is_empty() {
            continue;
        }

        let defline = if write_segments {
            projected_defline(transcript, &transcript.cds)
        } else {
            transcript.id.clone()
        };

        write_fasta_record(out, &defline, protein.as_bytes(), use_star_stop)
            .map_err(|err| CompatError::new(format!("Error writing FASTA: {err}\n"), 1))?;
    }
    Ok(())
}

pub fn spliced_sequence(
    transcript: &Transcript,
    segments: &[Segment],
    genome: &Genome,
) -> Result<Vec<u8>, CompatError> {
    spliced_sequence_with_padding(transcript, segments, genome, 0)
}

pub fn spliced_sequence_with_padding(
    transcript: &Transcript,
    segments: &[Segment],
    genome: &Genome,
    padding: u64,
) -> Result<Vec<u8>, CompatError> {
    let chrom = genome.get(&transcript.seqid).ok_or_else(|| {
        CompatError::new(
            format!(
                "Error: couldn't find genomic sequence {}\n",
                transcript.seqid
            ),
            1,
        )
    })?;

    let mut seq = Vec::new();
    let mut padded_segments = segments.iter().cloned().collect::<Vec<_>>();
    if padding > 0 && !padded_segments.is_empty() {
        let last_index = padded_segments.len() - 1;
        padded_segments[0].start = padded_segments[0].start.saturating_sub(padding);
        padded_segments[last_index].end =
            (padded_segments[last_index].end + padding).min(chrom.len() as u64);
    }

    for segment in &padded_segments {
        let start = segment.start.saturating_sub(1) as usize;
        let end = segment.end as usize;
        if start >= end || end > chrom.len() {
            return Err(CompatError::new(
                "Error: genomic segment outside sequence bounds\n",
                1,
            ));
        }
        seq.extend_from_slice(&chrom[start..end]);
    }

    if transcript.strand == '-' {
        reverse_complement(&mut seq);
    }

    Ok(seq)
}

pub fn unspliced_sequence(
    transcript: &Transcript,
    genome: &Genome,
    padding: u64,
) -> Result<Vec<u8>, CompatError> {
    let chrom = genome.get(&transcript.seqid).ok_or_else(|| {
        CompatError::new(
            format!(
                "Error: couldn't find genomic sequence {}\n",
                transcript.seqid
            ),
            1,
        )
    })?;

    let start = transcript.start.saturating_sub(padding).max(1);
    let end = (transcript.end + padding).min(chrom.len() as u64);
    let mut seq = chrom[(start - 1) as usize..end as usize].to_vec();
    if transcript.strand == '-' {
        reverse_complement(&mut seq);
    }
    Ok(seq)
}

pub fn reverse_complement(seq: &mut Vec<u8>) {
    seq.reverse();
    for base in seq {
        *base = match *base {
            b'A' | b'a' => b'T',
            b'C' | b'c' => b'G',
            b'G' | b'g' => b'C',
            b'T' | b't' => b'A',
            b'U' | b'u' => b'A',
            other => other.to_ascii_uppercase(),
        };
    }
}

pub fn translate(cds: &[u8]) -> String {
    let mut protein = cds
        .chunks(3)
        .filter(|codon| codon.len() == 3)
        .map(translate_codon)
        .collect::<String>();

    if protein.ends_with('.') {
        protein.pop();
    }

    protein
}

fn transcript_defline(
    transcript: &Transcript,
    genome: &Genome,
    padding: u64,
    suppress_cds: bool,
    write_segments: bool,
) -> Result<String, CompatError> {
    let mut defline = transcript.id.clone();

    if !suppress_cds {
        if let Some((start, end)) = projected_cds_span(transcript) {
            let prefix_padding = transcript_prefix_padding(transcript, genome, padding)?;
            defline.push_str(&format!(
                " CDS={}-{}",
                start + prefix_padding,
                end + prefix_padding
            ));
        }
    }

    if write_segments {
        defline.push_str(&spliced_segments_suffix(transcript, padding, genome)?);
    }

    Ok(defline)
}

fn spliced_segments_suffix(
    transcript: &Transcript,
    padding: u64,
    genome: &Genome,
) -> Result<String, CompatError> {
    let chrom = genome.get(&transcript.seqid).ok_or_else(|| {
        CompatError::new(
            format!(
                "Error: couldn't find genomic sequence {}\n",
                transcript.seqid
            ),
            1,
        )
    })?;

    let chrom_len = chrom.len() as u64;
    let left_padding = transcript.start.saturating_sub(1).min(padding);
    let right_padding = chrom_len.saturating_sub(transcript.end).min(padding);

    let mut suffix = format!(
        " loc:{}|{}-{}|{} exons:{}",
        transcript.seqid,
        transcript.start,
        transcript.end,
        transcript.strand,
        transcript.exon_list()
    );

    if padding > 0 {
        suffix.push_str(&format!(" padding:{}|{}", left_padding, right_padding));
    }

    let mut offset = 1u64;
    let mut parts = Vec::new();
    for (index, exon) in transcript.exons.iter().enumerate() {
        let extra_left = if index == 0 { left_padding } else { 0 };
        let extra_right = if index + 1 == transcript.exons.len() {
            right_padding
        } else {
            0
        };
        let len = segment_len(exon) + extra_left + extra_right;
        parts.push(format!("{offset}-{}", offset + len - 1));
        offset += len;
    }

    suffix.push_str(&format!(" segs:{}", parts.join(",")));
    Ok(suffix)
}

fn projected_defline(transcript: &Transcript, segments: &[Segment]) -> String {
    let mut offset = 1u64;
    let parts = projected_cds_segments(transcript, segments)
        .unwrap_or_default()
        .into_iter()
        .map(|segment| {
            let len = segment_len(&segment);
            let part = format!("{offset}-{}", offset + len - 1);
            offset += len;
            part
        })
        .collect::<Vec<_>>()
        .join(",");

    format!(
        "{} loc:{}({}){}-{} segs:{}",
        transcript.id, transcript.seqid, transcript.strand, transcript.start, transcript.end, parts
    )
}

fn transcript_prefix_padding(
    transcript: &Transcript,
    genome: &Genome,
    padding: u64,
) -> Result<u64, CompatError> {
    if padding == 0 {
        return Ok(0);
    }

    let chrom = genome.get(&transcript.seqid).ok_or_else(|| {
        CompatError::new(
            format!(
                "Error: couldn't find genomic sequence {}\n",
                transcript.seqid
            ),
            1,
        )
    })?;

    let chrom_len = chrom.len() as u64;
    let prefix = if transcript.strand == '-' {
        chrom_len.saturating_sub(transcript.end)
    } else {
        transcript.start.saturating_sub(1)
    };

    Ok(prefix.min(padding))
}

pub fn write_fasta_record<W: Write>(
    out: &mut W,
    defline: &str,
    seq: &[u8],
    use_star_stop: bool,
) -> io::Result<()> {
    writeln!(out, ">{defline}")?;
    for chunk in seq.chunks(70) {
        if use_star_stop {
            let converted = chunk
                .iter()
                .map(|base| if *base == b'.' { b'*' } else { *base })
                .collect::<Vec<_>>();
            out.write_all(&converted)?;
        } else {
            out.write_all(chunk)?;
        }
        out.write_all(b"\n")?;
    }
    Ok(())
}

pub fn translate_codon(codon: &[u8]) -> char {
    match upper_codon(codon).as_str() {
        "TTT" | "TTC" => 'F',
        "TTA" | "TTG" | "CTT" | "CTC" | "CTA" | "CTG" => 'L',
        "ATT" | "ATC" | "ATA" => 'I',
        "ATG" => 'M',
        "GTT" | "GTC" | "GTA" | "GTG" => 'V',
        "TCT" | "TCC" | "TCA" | "TCG" | "AGT" | "AGC" => 'S',
        "CCT" | "CCC" | "CCA" | "CCG" => 'P',
        "ACT" | "ACC" | "ACA" | "ACG" => 'T',
        "GCT" | "GCC" | "GCA" | "GCG" => 'A',
        "TAT" | "TAC" => 'Y',
        "TAA" | "TAG" | "TGA" => '.',
        "CAT" | "CAC" => 'H',
        "CAA" | "CAG" => 'Q',
        "AAT" | "AAC" => 'N',
        "AAA" | "AAG" => 'K',
        "GAT" | "GAC" => 'D',
        "GAA" | "GAG" => 'E',
        "TGT" | "TGC" => 'C',
        "TGG" => 'W',
        "CGT" | "CGC" | "CGA" | "CGG" | "AGA" | "AGG" => 'R',
        "GGT" | "GGC" | "GGA" | "GGG" => 'G',
        _ => 'X',
    }
}

fn upper_codon(codon: &[u8]) -> String {
    codon
        .iter()
        .map(|base| (*base as char).to_ascii_uppercase())
        .collect()
}

fn finalize_record(
    genome: &mut Genome,
    index_entries: &mut Vec<FastaIndexEntry>,
    current_id: &mut Option<String>,
    current_seq: &mut Vec<u8>,
    current_offset: usize,
    current_line_bases: usize,
    current_line_width: usize,
) {
    let Some(id) = current_id.take() else {
        return;
    };

    let len = current_seq.len();
    genome.insert(id.clone(), std::mem::take(current_seq));
    index_entries.push(FastaIndexEntry {
        id,
        len,
        offset: current_offset,
        line_bases: current_line_bases,
        line_width: current_line_width.max(current_line_bases),
    });
}

fn trim_line_ending(line: &[u8]) -> &[u8] {
    let line = line.strip_suffix(b"\n").unwrap_or(line);
    line.strip_suffix(b"\r").unwrap_or(line)
}

fn parse_record_id(header: &[u8]) -> String {
    String::from_utf8_lossy(header)
        .split_whitespace()
        .next()
        .unwrap_or("")
        .to_owned()
}

fn write_fasta_index(path: &Path, entries: &[FastaIndexEntry]) -> Result<(), CompatError> {
    let fai_path = fasta_index_path(path);
    let existed = fai_path.exists();
    let mut text = String::new();
    for entry in entries {
        let _ = writeln!(
            &mut text,
            "{}\t{}\t{}\t{}\t{}",
            entry.id, entry.len, entry.offset, entry.line_bases, entry.line_width
        );
    }

    fs::write(&fai_path, text).map_err(|_| {
        CompatError::new(
            format!(
                "Error: couldn't write FASTA index file {}\n",
                fai_path.display()
            ),
            1,
        )
    })?;

    if !existed {
        eprint!("FASTA index file {} created.\n", fai_path.display());
    }

    Ok(())
}

fn fasta_index_path(path: &Path) -> PathBuf {
    let mut file_name = path.as_os_str().to_os_string();
    file_name.push(".fai");
    PathBuf::from(file_name)
}

fn projected_cds_span(transcript: &Transcript) -> Option<(u64, u64)> {
    if transcript.cds.is_empty() || transcript.exons.is_empty() {
        return None;
    }

    let (first_index, first_segment) = transcript_first_cds_segment(transcript)?;
    let (coding_start, coding_end) =
        trimmed_segment_bounds(transcript, first_segment, first_index)?;
    let coding_pos = if transcript.strand == '-' {
        coding_end
    } else {
        coding_start
    };

    let mut transcript_offset = 1u64;
    for exon in transcript_ordered_segments(transcript, &transcript.exons) {
        if coding_pos >= exon.start && coding_pos <= exon.end {
            let start = if transcript.strand == '-' {
                transcript_offset + (exon.end - coding_pos)
            } else {
                transcript_offset + (coding_pos - exon.start)
            };
            let cds_len = cds_sequence_len(transcript)?;
            return Some((start, start + cds_len - 1));
        }
        transcript_offset += segment_len(exon);
    }

    None
}

fn transcript_ordered_segments<'a>(
    transcript: &Transcript,
    segments: &'a [Segment],
) -> Vec<&'a Segment> {
    let mut ordered = segments.iter().collect::<Vec<_>>();
    if transcript.strand == '-' {
        ordered.reverse();
    }
    ordered
}

fn segment_len(segment: &Segment) -> u64 {
    segment.end - segment.start + 1
}

fn spliced_cds_sequence(transcript: &Transcript, genome: &Genome) -> Result<Vec<u8>, CompatError> {
    let chrom = genome.get(&transcript.seqid).ok_or_else(|| {
        CompatError::new(
            format!(
                "Error: couldn't find genomic sequence {}\n",
                transcript.seqid
            ),
            1,
        )
    })?;

    let mut seq = Vec::new();
    for (index, segment) in transcript.cds.iter().enumerate() {
        let (trimmed_start, trimmed_end) = trimmed_segment_bounds(transcript, segment, index)
            .ok_or_else(|| CompatError::new("Error: invalid CDS phase trimming\n", 1))?;
        let start = trimmed_start.saturating_sub(1) as usize;
        let end = trimmed_end as usize;
        if start >= end || end > chrom.len() {
            return Err(CompatError::new(
                "Error: genomic segment outside sequence bounds\n",
                1,
            ));
        }
        seq.extend_from_slice(&chrom[start..end]);
    }

    if transcript.strand == '-' {
        reverse_complement(&mut seq);
    }

    Ok(seq)
}

fn trimmed_segment_bounds(
    transcript: &Transcript,
    segment: &Segment,
    index: usize,
) -> Option<(u64, u64)> {
    let trim_phase = if transcript.strand == '-' {
        index + 1 == transcript.cds.len()
    } else {
        index == 0
    };

    if trim_phase {
        let phase = segment.phase.parse::<u64>().ok()?;
        if transcript.strand == '-' {
            let end = segment.end.checked_sub(phase)?;
            (segment.start <= end).then_some((segment.start, end))
        } else {
            let start = segment.start.checked_add(phase)?;
            (start <= segment.end).then_some((start, segment.end))
        }
    } else {
        Some((segment.start, segment.end))
    }
}

fn transcript_first_cds_segment(transcript: &Transcript) -> Option<(usize, &Segment)> {
    if transcript.strand == '-' {
        transcript.cds.iter().enumerate().next_back()
    } else {
        transcript.cds.iter().enumerate().next()
    }
}

fn cds_sequence_len(transcript: &Transcript) -> Option<u64> {
    transcript
        .cds
        .iter()
        .enumerate()
        .map(|(index, segment)| {
            trimmed_segment_bounds(transcript, segment, index).map(|(start, end)| end - start + 1)
        })
        .sum()
}

fn projected_cds_segments<'a>(
    transcript: &Transcript,
    segments: &'a [Segment],
) -> Option<Vec<Segment>> {
    let mut projected = segments
        .iter()
        .enumerate()
        .map(|(index, segment)| {
            let (start, end) = trimmed_segment_bounds(transcript, segment, index)?;
            Some(Segment {
                start,
                end,
                phase: segment.phase.clone(),
            })
        })
        .collect::<Option<Vec<_>>>()?;

    if transcript.strand == '-' {
        projected.reverse();
    }

    Some(projected)
}

struct FastaIndexEntry {
    id: String,
    len: usize,
    offset: usize,
    line_bases: usize,
    line_width: usize,
}
