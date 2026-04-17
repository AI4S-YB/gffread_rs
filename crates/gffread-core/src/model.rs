use std::collections::BTreeMap;

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Annotation {
    pub transcripts: Vec<Transcript>,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Locus {
    pub seqid: String,
    pub start: u64,
    pub end: u64,
    pub strand: char,
    pub id: String,
    pub genes: Vec<String>,
    pub gene_ids: Vec<String>,
    pub transcripts: Vec<String>,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Transcript {
    pub seqid: String,
    pub source: String,
    pub feature: String,
    pub start: u64,
    pub end: u64,
    pub strand: char,
    pub id: String,
    pub gene_id: Option<String>,
    pub gene_name: Option<String>,
    pub attrs: BTreeMap<String, String>,
    pub locus: Option<String>,
    pub exons: Vec<Segment>,
    pub cds: Vec<Segment>,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Segment {
    pub start: u64,
    pub end: u64,
    pub phase: String,
}

impl Transcript {
    pub fn exon_list(&self) -> String {
        self.exons
            .iter()
            .map(|segment| format!("{}-{}", segment.start, segment.end))
            .collect::<Vec<_>>()
            .join(",")
    }

    pub fn covlen(&self) -> u64 {
        self.exons
            .iter()
            .map(|segment| segment.end - segment.start + 1)
            .sum()
    }

    pub fn max_intron_len(&self) -> u64 {
        self.exons
            .windows(2)
            .map(|segments| segments[1].start.saturating_sub(segments[0].end + 1))
            .max()
            .unwrap_or(0)
    }

    pub fn is_pseudo(&self) -> bool {
        self.attrs.iter().any(|(key, value)| {
            key.eq_ignore_ascii_case("pseudo") || value.to_ascii_lowercase().contains("pseudo")
        })
    }
}
