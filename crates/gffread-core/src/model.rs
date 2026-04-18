#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Annotation {
    pub transcripts: Vec<Transcript>,
    pub genes: Vec<Gene>,
    pub ref_order: Vec<String>,
    pub header_comments: Vec<String>,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Attr {
    pub key: String,
    pub value: String,
}

#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct Attrs {
    items: Vec<Attr>,
}

impl Attrs {
    pub fn new() -> Self {
        Self { items: Vec::new() }
    }

    pub fn push_unique(&mut self, key: String, value: String) {
        if self.items.iter().any(|item| item.key == key) {
            return;
        }
        self.items.push(Attr { key, value });
    }

    pub fn insert_or_replace(&mut self, key: String, value: String) {
        if let Some(item) = self.items.iter_mut().find(|item| item.key == key) {
            item.value = value;
            return;
        }
        self.items.push(Attr { key, value });
    }

    pub fn get(&self, key: &str) -> Option<&str> {
        self.items
            .iter()
            .find(|item| item.key == key && !item.value.is_empty())
            .map(|item| item.value.as_str())
    }

    pub fn contains_key(&self, key: &str) -> bool {
        self.items.iter().any(|item| item.key == key)
    }

    pub fn iter(&self) -> impl Iterator<Item = &Attr> {
        self.items.iter()
    }

    pub fn is_empty(&self) -> bool {
        self.items.is_empty()
    }

    pub fn remove_all(&mut self) {
        self.items.clear();
    }
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
    pub order: usize,
    pub gene_id: Option<String>,
    pub gene_name: Option<String>,
    pub attrs: Attrs,
    pub locus: Option<String>,
    pub exons: Vec<Segment>,
    pub cds: Vec<Segment>,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Gene {
    pub seqid: String,
    pub source: String,
    pub feature: String,
    pub start: u64,
    pub end: u64,
    pub strand: char,
    pub id: String,
    pub order: usize,
    pub attrs: Attrs,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Segment {
    pub start: u64,
    pub end: u64,
    pub score: String,
    pub phase: String,
    pub attrs: Attrs,
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
        self.attrs.iter().any(|attr| {
            attr.key.eq_ignore_ascii_case("pseudo")
                || attr.value.to_ascii_lowercase().contains("pseudo")
        })
    }

    pub fn has_cds(&self) -> bool {
        !self.cds.is_empty()
    }

    pub fn cds_phase(&self) -> Option<char> {
        self.cds
            .first()
            .and_then(|segment| segment.phase.chars().next())
            .filter(|phase| matches!(phase, '0' | '1' | '2'))
    }
}

impl Gene {
    pub fn gene_name(&self) -> Option<&str> {
        self.attrs
            .get("gene_name")
            .or_else(|| self.attrs.get("gene_sym"))
            .or_else(|| self.attrs.get("gene"))
            .or_else(|| self.attrs.get("genesymbol"))
    }
}
