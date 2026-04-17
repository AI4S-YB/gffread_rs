use std::collections::BTreeMap;

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Annotation {
    pub transcripts: Vec<Transcript>,
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
}
