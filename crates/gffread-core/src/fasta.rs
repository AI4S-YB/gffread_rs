use std::cell::RefCell;
use std::collections::BTreeMap;
use std::fmt::Write as _;
use std::fs::{self, File};
use std::io::{self, BufRead, BufReader, Read, Seek, SeekFrom, Write};
use std::path::{Path, PathBuf};
use std::rc::Rc;

use crate::compat::CompatError;
use crate::model::{Annotation, Attrs, Segment, Transcript};

pub struct Genome {
    index: BTreeMap<String, FastaIndexEntry>,
    file: RefCell<File>,
    cache: RefCell<BTreeMap<String, Rc<Vec<u8>>>>,
    window: RefCell<Option<SequenceWindow>>,
}

impl Genome {
    pub fn get(&self, seqid: &str) -> Option<Rc<Vec<u8>>> {
        if let Some(seq) = self.cache.borrow().get(seqid).cloned() {
            return Some(seq);
        }

        let len = self.len_of(seqid)?;
        let seq = Rc::new(self.read_range(seqid, 1, len).ok()?);
        self.cache
            .borrow_mut()
            .insert(seqid.to_owned(), seq.clone());
        Some(seq)
    }

    pub fn len_of(&self, seqid: &str) -> Option<u64> {
        self.index.get(seqid).map(|entry| entry.len as u64)
    }

    pub fn read_range(&self, seqid: &str, start: u64, end: u64) -> Result<Vec<u8>, CompatError> {
        self.with_range(seqid, start, end, |seq| Ok(seq.to_vec()))
    }

    fn with_range<T, F>(&self, seqid: &str, start: u64, end: u64, f: F) -> Result<T, CompatError>
    where
        F: FnOnce(&[u8]) -> Result<T, CompatError>,
    {
        let entry = self
            .index
            .get(seqid)
            .ok_or_else(|| missing_sequence_error(seqid))?;
        if start == 0 || start > end || end > entry.len as u64 {
            return Err(CompatError::new(
                "Error: genomic segment outside sequence bounds\n",
                1,
            ));
        }

        let mut window_ref = self.window.borrow_mut();
        let cache_hit = window_ref.as_ref().is_some_and(|window| {
            window.seqid == seqid && start >= window.start && end <= window.end()
        });
        if !cache_hit {
            let mut file = self.file.borrow_mut();
            *window_ref = Some(load_window(seqid, entry, start, end, &mut file)?);
        }

        let window = window_ref
            .as_ref()
            .expect("window should be loaded before slicing");
        f(window.slice(start, end))
    }
}

pub fn load_genome(path: &Path) -> Result<Genome, CompatError> {
    let file = File::open(path).map_err(|_| {
        CompatError::new(
            format!(
                "Error: couldn't open genomic sequences file {}\n",
                path.display()
            ),
            1,
        )
    })?;

    let index = match read_fasta_index(path) {
        Ok(Some(index)) => index,
        _ => {
            let entries = build_fasta_index(path)?;
            write_fasta_index(path, &entries)?;
            entries
                .into_iter()
                .map(|entry| (entry.id.clone(), entry))
                .collect()
        }
    };

    Ok(Genome {
        index,
        file: RefCell::new(file),
        cache: RefCell::new(BTreeMap::new()),
        window: RefCell::new(None),
    })
}

fn missing_sequence_error(seqid: &str) -> CompatError {
    CompatError::new(
        format!("Error: couldn't find genomic sequence {seqid}\n"),
        1,
    )
}

fn chrom_len(genome: &Genome, seqid: &str) -> Result<u64, CompatError> {
    genome
        .len_of(seqid)
        .ok_or_else(|| missing_sequence_error(seqid))
}

#[derive(Clone)]
struct SequenceWindow {
    seqid: String,
    start: u64,
    seq: Vec<u8>,
}

impl SequenceWindow {
    fn end(&self) -> u64 {
        self.start + self.seq.len() as u64 - 1
    }

    fn slice(&self, start: u64, end: u64) -> &[u8] {
        let from = (start - self.start) as usize;
        let to = (end - self.start + 1) as usize;
        &self.seq[from..to]
    }
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
        let transcript = clipped_transcript(transcript, genome)?;
        let seq = spliced_sequence_with_padding(&transcript, &transcript.exons, genome, padding)?;
        let defline =
            transcript_defline(&transcript, genome, padding, suppress_cds, write_segments)?;
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
        let transcript = clipped_transcript(transcript, genome)?;
        let seq = unspliced_sequence(&transcript, genome, padding)?;
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
        let transcript = clipped_transcript(transcript, genome)?;
        if transcript.cds.is_empty() {
            continue;
        }

        let seq = spliced_cds_sequence(&transcript, genome)?;
        warn_short_cds(&transcript, &seq);
        let defline = if write_segments {
            projected_defline(&transcript, &transcript.cds)
        } else {
            transcript.id.clone()
        };
        write_fasta_record(out, &defline, &seq, false)
            .map_err(|err| CompatError::new(format!("Error writing FASTA: {err}\n"), 1))?;
    }
    Ok(())
}

fn warn_short_cds(transcript: &Transcript, seq: &[u8]) {
    if seq.len() >= 4 {
        return;
    }
    let Some(first_cds) = transcript.cds.first() else {
        return;
    };
    let Some(last_cds) = transcript.cds.last() else {
        return;
    };
    eprintln!(
        "Warning: CDS {}-{} too short for {}, check your data.",
        first_cds.start, last_cds.end, transcript.id
    );
}

pub fn write_protein_fasta<W: Write>(
    out: &mut W,
    annotation: &Annotation,
    genome: &Genome,
    use_star_stop: bool,
    write_segments: bool,
    warn_for_short_cds: bool,
) -> Result<(), CompatError> {
    for transcript in &annotation.transcripts {
        let transcript = clipped_transcript(transcript, genome)?;
        if transcript.cds.is_empty() {
            continue;
        }

        let cds = spliced_cds_sequence(&transcript, genome)?;
        if warn_for_short_cds {
            warn_short_cds(&transcript, &cds);
        }
        let protein = translate(&cds);
        if protein.is_empty() {
            continue;
        }

        let defline = if write_segments {
            projected_defline(&transcript, &transcript.cds)
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
    let chrom_len = chrom_len(genome, &transcript.seqid)?;
    let mut padded_segments = segments.to_vec();
    if padding > 0 && !padded_segments.is_empty() {
        let last_index = padded_segments.len() - 1;
        padded_segments[0].start = padded_segments[0].start.saturating_sub(padding);
        padded_segments[last_index].end =
            (padded_segments[last_index].end + padding).min(chrom_len);
    }

    if padded_segments.is_empty() {
        return Ok(Vec::new());
    }

    let span_start = padded_segments[0].start;
    let span_end = padded_segments[padded_segments.len() - 1]
        .end
        .min(chrom_len);
    if span_start == 0 || span_start > span_end {
        return Err(CompatError::new(
            "Error: genomic segment outside sequence bounds\n",
            1,
        ));
    }

    genome.with_range(&transcript.seqid, span_start, span_end, |span| {
        let mut seq = Vec::with_capacity(
            padded_segments
                .iter()
                .map(|segment| segment.end.min(chrom_len).saturating_sub(segment.start) + 1)
                .sum::<u64>() as usize,
        );

        for segment in &padded_segments {
            let clipped_end = segment.end.min(chrom_len);
            if segment.start == 0 || segment.start > clipped_end {
                return Err(CompatError::new(
                    "Error: genomic segment outside sequence bounds\n",
                    1,
                ));
            }
            let rel_start = (segment.start - span_start) as usize;
            let rel_end = (clipped_end - span_start + 1) as usize;
            seq.extend_from_slice(&span[rel_start..rel_end]);
        }

        if transcript.strand == '-' {
            reverse_complement(&mut seq);
        }

        Ok(seq)
    })
}

pub fn unspliced_sequence(
    transcript: &Transcript,
    genome: &Genome,
    padding: u64,
) -> Result<Vec<u8>, CompatError> {
    let chrom_len = chrom_len(genome, &transcript.seqid)?;
    let start = transcript.start.saturating_sub(padding).max(1);
    let end = (transcript.end + padding).min(chrom_len);
    genome.with_range(&transcript.seqid, start, end, |span| {
        let mut seq = span.to_vec();
        if transcript.strand == '-' {
            reverse_complement(&mut seq);
        }
        Ok(seq)
    })
}

pub fn reverse_complement(seq: &mut Vec<u8>) {
    seq.reverse();
    for base in seq {
        *base = match *base {
            b'A' => b'T',
            b'a' => b't',
            b'C' => b'G',
            b'c' => b'g',
            b'G' => b'C',
            b'g' => b'c',
            b'T' => b'A',
            b't' => b'a',
            b'U' => b'A',
            b'u' => b'a',
            other => other,
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
    let chrom_len = chrom_len(genome, &transcript.seqid)?;
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

    let chrom_len = chrom_len(genome, &transcript.seqid)?;
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
    let possibilities = codon_possibilities(codon);
    if possibilities.is_empty() {
        return 'X';
    }

    let mut aa_options = Vec::new();
    for candidate in possibilities {
        let aa = translate_resolved_codon(&candidate);
        if !aa_options.contains(&aa) {
            aa_options.push(aa);
        }
    }

    if aa_options.len() == 1 {
        return aa_options[0];
    }
    aa_options.sort_unstable();
    match aa_options.as_slice() {
        ['D', 'N'] => 'B',
        ['E', 'Q'] => 'Z',
        _ => 'X',
    }
}

fn translate_resolved_codon(codon: &[u8; 3]) -> char {
    match std::str::from_utf8(codon).unwrap_or("") {
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

fn codon_possibilities(codon: &[u8]) -> Vec<[u8; 3]> {
    if codon.len() != 3 {
        return Vec::new();
    }

    let mut result = Vec::new();
    for first in possible_bases(codon[0]) {
        for second in possible_bases(codon[1]) {
            for third in possible_bases(codon[2]) {
                result.push([*first, *second, *third]);
            }
        }
    }
    result
}

fn possible_bases(base: u8) -> &'static [u8] {
    match base.to_ascii_uppercase() {
        b'A' => b"A",
        b'C' => b"C",
        b'G' => b"G",
        b'T' | b'U' => b"T",
        b'R' => b"AG",
        b'Y' => b"CT",
        b'S' => b"CG",
        b'W' => b"AT",
        b'K' => b"GT",
        b'M' => b"AC",
        b'B' => b"CGT",
        b'D' => b"AGT",
        b'H' => b"ACT",
        b'V' => b"ACG",
        b'N' | b'X' => b"ACGT",
        _ => b"",
    }
}

fn clipped_transcript(transcript: &Transcript, genome: &Genome) -> Result<Transcript, CompatError> {
    let chrom_len = chrom_len(genome, &transcript.seqid)?;
    let mut clipped = transcript.clone();
    clip_segments_to_chrom(&mut clipped.exons, chrom_len)?;
    clip_segments_to_chrom(&mut clipped.cds, chrom_len)?;
    if let Some(first_exon) = clipped.exons.first() {
        clipped.start = first_exon.start;
    }
    if let Some(last_exon) = clipped.exons.last() {
        clipped.end = last_exon.end;
    }
    Ok(clipped)
}

fn clip_segments_to_chrom(segments: &mut [Segment], chrom_len: u64) -> Result<(), CompatError> {
    for segment in segments {
        if segment.end > chrom_len {
            segment.end = chrom_len;
        }
        if segment.start == 0 || segment.start > segment.end {
            return Err(CompatError::new(
                "Error: genomic segment outside sequence bounds\n",
                1,
            ));
        }
    }
    Ok(())
}

fn finalize_index_entry(
    index_entries: &mut Vec<FastaIndexEntry>,
    current_id: &mut Option<String>,
    current_len: &mut usize,
    current_offset: usize,
    current_line_bases: usize,
    current_line_width: usize,
) {
    let Some(id) = current_id.take() else {
        return;
    };

    index_entries.push(FastaIndexEntry {
        id,
        len: *current_len,
        offset: current_offset,
        line_bases: current_line_bases,
        line_width: current_line_width.max(current_line_bases),
    });
    *current_len = 0;
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

fn read_fasta_index(path: &Path) -> Result<Option<BTreeMap<String, FastaIndexEntry>>, CompatError> {
    let fai_path = fasta_index_path(path);
    let file = match File::open(&fai_path) {
        Ok(file) => file,
        Err(err) if err.kind() == io::ErrorKind::NotFound => return Ok(None),
        Err(_) => return Ok(None),
    };

    let mut index = BTreeMap::new();
    for line in BufReader::new(file).lines() {
        let Ok(line) = line else {
            return Ok(None);
        };
        if line.trim().is_empty() {
            continue;
        }
        let fields = line.split('\t').collect::<Vec<_>>();
        if fields.len() != 5 {
            return Ok(None);
        }
        let (Ok(len), Ok(offset), Ok(line_bases), Ok(line_width)) = (
            fields[1].parse::<usize>(),
            fields[2].parse::<usize>(),
            fields[3].parse::<usize>(),
            fields[4].parse::<usize>(),
        ) else {
            return Ok(None);
        };
        index.insert(
            fields[0].to_owned(),
            FastaIndexEntry {
                id: fields[0].to_owned(),
                len,
                offset,
                line_bases,
                line_width,
            },
        );
    }

    if index.is_empty() {
        Ok(None)
    } else {
        Ok(Some(index))
    }
}

fn build_fasta_index(path: &Path) -> Result<Vec<FastaIndexEntry>, CompatError> {
    let file = File::open(path).map_err(|_| {
        CompatError::new(
            format!(
                "Error: couldn't open genomic sequences file {}\n",
                path.display()
            ),
            1,
        )
    })?;

    let mut index_entries = Vec::new();
    let mut current_id = None::<String>;
    let mut current_len = 0usize;
    let mut current_offset = 0usize;
    let mut current_line_bases = 0usize;
    let mut current_line_width = 0usize;
    let mut position = 0usize;
    let mut reader = BufReader::with_capacity(1024 * 1024, file);
    let mut raw_line = Vec::new();

    loop {
        raw_line.clear();
        let bytes_read = reader
            .read_until(b'\n', &mut raw_line)
            .map_err(|_| CompatError::new("Error: couldn't read genomic sequence file\n", 1))?;
        if bytes_read == 0 {
            break;
        }

        let line = trim_line_ending(&raw_line);
        if let Some(rest) = line.strip_prefix(b">") {
            finalize_index_entry(
                &mut index_entries,
                &mut current_id,
                &mut current_len,
                current_offset,
                current_line_bases,
                current_line_width,
            );

            current_id = Some(parse_record_id(rest));
            current_offset = position + raw_line.len();
            current_line_bases = 0;
            current_line_width = 0;
        } else if current_id.is_some() {
            let bases_len = line
                .iter()
                .copied()
                .filter(|byte| !byte.is_ascii_whitespace())
                .count();

            if bases_len > 0 {
                if current_line_bases == 0 {
                    current_line_bases = bases_len;
                    current_line_width = raw_line.len();
                }
                current_len += bases_len;
            }
        }

        position += raw_line.len();
    }

    finalize_index_entry(
        &mut index_entries,
        &mut current_id,
        &mut current_len,
        current_offset,
        current_line_bases,
        current_line_width,
    );

    Ok(index_entries)
}

fn read_range_from_index(
    file: &mut File,
    entry: &FastaIndexEntry,
    start: u64,
    end: u64,
) -> Result<Vec<u8>, CompatError> {
    if start == 0 || start > end || end > entry.len as u64 {
        return Err(CompatError::new(
            "Error: genomic segment outside sequence bounds\n",
            1,
        ));
    }
    if entry.line_bases == 0 {
        return Err(CompatError::new(
            "Error: FASTA index contains an empty sequence line width\n",
            1,
        ));
    }

    let start0 = start - 1;
    let end0 = end - 1;
    let start_line = start0 / entry.line_bases as u64;
    let start_col = start0 % entry.line_bases as u64;
    let end_line = end0 / entry.line_bases as u64;
    let end_col = end0 % entry.line_bases as u64;
    let file_start = entry.offset as u64 + start_line * entry.line_width as u64 + start_col;
    let file_end = entry.offset as u64 + end_line * entry.line_width as u64 + end_col + 1;
    let span_len = (file_end - file_start) as usize;

    file.seek(SeekFrom::Start(file_start))
        .map_err(|_| CompatError::new("Error: couldn't seek genomic sequence file\n", 1))?;

    let mut raw = vec![0u8; span_len];
    file.read_exact(&mut raw)
        .map_err(|_| CompatError::new("Error: couldn't read genomic sequence file\n", 1))?;

    let expected_len = (end - start + 1) as usize;
    let mut seq = Vec::with_capacity(expected_len);
    seq.extend(raw.into_iter().filter(|byte| !byte.is_ascii_whitespace()));

    if seq.len() != expected_len {
        return Err(CompatError::new(
            "Error: genomic sequence truncated while reading FASTA\n",
            1,
        ));
    }

    Ok(seq)
}

fn load_window(
    seqid: &str,
    entry: &FastaIndexEntry,
    start: u64,
    end: u64,
    file: &mut File,
) -> Result<SequenceWindow, CompatError> {
    Ok(SequenceWindow {
        seqid: seqid.to_owned(),
        start,
        seq: read_range_from_index(file, entry, start, end)?,
    })
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
        eprintln!("FASTA index file {} created.", fai_path.display());
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

    let projected = projected_cds_segments(transcript, &transcript.cds)?;
    let mut offset = 1u64;
    let mut start = None;
    let mut end = None;
    for exon in transcript_ordered_segments(transcript, &transcript.exons) {
        let exon_start = offset;
        let exon_end = offset + segment_len(exon) - 1;
        for segment in &projected {
            let overlap_start = segment.start.max(exon.start);
            let overlap_end = segment.end.min(exon.end);
            if overlap_start > overlap_end {
                continue;
            }
            let seg_start = if transcript.strand == '-' {
                exon_start + (exon.end - overlap_end)
            } else {
                exon_start + (overlap_start - exon.start)
            };
            let seg_end = seg_start + (overlap_end - overlap_start);
            start.get_or_insert(seg_start);
            end = Some(seg_end);
        }
        offset = exon_end + 1;
    }

    start.zip(end)
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
    let chrom_len = chrom_len(genome, &transcript.seqid)?;
    let span_start = transcript.start;
    let span_end = transcript.end.min(chrom_len);
    if span_start == 0 || span_start > span_end {
        return Err(CompatError::new(
            "Error: genomic segment outside sequence bounds\n",
            1,
        ));
    }
    genome.with_range(&transcript.seqid, span_start, span_end, |span| {
        let mut seq = Vec::new();
        for (index, segment) in transcript.cds.iter().enumerate() {
            let (trimmed_start, trimmed_end) =
                trimmed_segment_bounds(transcript, segment, index)
                    .ok_or_else(|| CompatError::new("Error: invalid CDS phase trimming\n", 1))?;
            let trimmed_end = trimmed_end.min(chrom_len);
            if trimmed_start == 0 || trimmed_start > trimmed_end {
                return Err(CompatError::new(
                    "Error: genomic segment outside sequence bounds\n",
                    1,
                ));
            }
            let rel_start = (trimmed_start - span_start) as usize;
            let rel_end = (trimmed_end - span_start + 1) as usize;
            seq.extend_from_slice(&span[rel_start..rel_end]);
        }

        if transcript.strand == '-' {
            reverse_complement(&mut seq);
        }

        Ok(seq)
    })
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

fn projected_cds_segments(transcript: &Transcript, segments: &[Segment]) -> Option<Vec<Segment>> {
    let mut projected = segments
        .iter()
        .enumerate()
        .map(|(index, segment)| {
            let (start, end) = trimmed_segment_bounds(transcript, segment, index)?;
            Some(Segment {
                start,
                end,
                score: ".".to_owned(),
                phase: segment.phase.clone(),
                attrs: Attrs::new(),
            })
        })
        .collect::<Option<Vec<_>>>()?;

    if transcript.strand == '-' {
        projected.reverse();
    }

    Some(projected)
}

#[derive(Clone)]
struct FastaIndexEntry {
    id: String,
    len: usize,
    offset: usize,
    line_bases: usize,
    line_width: usize,
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::time::{SystemTime, UNIX_EPOCH};

    fn temp_path(name: &str) -> PathBuf {
        std::env::temp_dir().join(format!(
            "gffread-rs-{name}-{}-{}",
            std::process::id(),
            SystemTime::now()
                .duration_since(UNIX_EPOCH)
                .expect("clock should be after epoch")
                .as_nanos()
        ))
    }

    #[test]
    fn fasta_index_roundtrip_reads_expected_subsequence() {
        let fasta_path = temp_path("fasta-index.fa");
        let fai_path = PathBuf::from(format!("{}.fai", fasta_path.display()));
        fs::write(&fasta_path, b">chr1\nAaCcGgTt\nAaCcGgTt\n>chr2\nTTttGGgg\n")
            .expect("fixture FASTA should be written");

        let genome = load_genome(&fasta_path).expect("genome should load");
        let chr1 = genome.get("chr1").expect("chr1 should be present");
        let chr2 = genome.get("chr2").expect("chr2 should be present");

        assert_eq!(&*chr1, b"AaCcGgTtAaCcGgTt");
        assert_eq!(&*chr2, b"TTttGGgg");
        assert!(fai_path.exists(), ".fai should be created");

        let _ = fs::remove_file(&fasta_path);
        let _ = fs::remove_file(&fai_path);
    }

    #[test]
    fn fasta_index_file_is_reused_for_lazy_sequence_loading() {
        let fasta_path = temp_path("fasta-index-existing.fa");
        let fai_path = PathBuf::from(format!("{}.fai", fasta_path.display()));
        fs::write(&fasta_path, b">chr1\nAaCcGgTt\nAaCcGgTt\n>chr2\nTTttGGgg\n")
            .expect("fixture FASTA should be written");
        fs::write(&fai_path, b"chr1\t16\t6\t8\t9\nchr2\t8\t30\t8\t9\n")
            .expect("fixture FASTA index should be written");

        let genome = load_genome(&fasta_path).expect("genome should load from existing index");
        let chr2 = genome.get("chr2").expect("chr2 should be present");
        assert_eq!(&*chr2, b"TTttGGgg");

        let _ = fs::remove_file(&fasta_path);
        let _ = fs::remove_file(&fai_path);
    }

    #[test]
    fn fasta_index_reads_requested_interval_across_wrapped_lines() {
        let fasta_path = temp_path("fasta-index-interval.fa");
        let fai_path = PathBuf::from(format!("{}.fai", fasta_path.display()));
        fs::write(&fasta_path, b">chr1\nACGTacgt\nTTGGccaa\n")
            .expect("fixture FASTA should be written");

        let genome = load_genome(&fasta_path).expect("genome should load");
        let interval = genome
            .read_range("chr1", 7, 12)
            .expect("interval should be readable");

        assert_eq!(interval, b"gtTTGG");

        let _ = fs::remove_file(&fasta_path);
        let _ = fs::remove_file(&fai_path);
    }

    #[test]
    fn translate_codon_matches_gclib_ambiguous_codon_table() {
        assert_eq!(translate_codon(b"GGN"), 'G');
        assert_eq!(translate_codon(b"ACN"), 'T');
        assert_eq!(translate_codon(b"CGN"), 'R');
        assert_eq!(translate_codon(b"TAR"), '.');
        assert_eq!(translate_codon(b"RAY"), 'B');
        assert_eq!(translate_codon(b"SAR"), 'Z');
        assert_eq!(translate_codon(b"NNN"), 'X');
    }
}
