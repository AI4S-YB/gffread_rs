use std::collections::BTreeSet;
use std::fs;
use std::path::Path;

use crate::compat::CompatError;
use crate::model::{Annotation, Transcript};
use crate::options::{IdFilterMode, RuntimeOptions};

pub fn apply_filters(
    annotation: &Annotation,
    options: &RuntimeOptions,
) -> Result<Annotation, CompatError> {
    let id_set = match &options.id_filter {
        Some(filter) => Some(load_ids(&filter.path)?),
        None => None,
    };

    let transcripts = annotation
        .transcripts
        .iter()
        .filter(|transcript| passes_filters(transcript, options, id_set.as_ref()))
        .cloned()
        .collect();

    Ok(Annotation { transcripts })
}

fn passes_filters(
    transcript: &Transcript,
    options: &RuntimeOptions,
    id_set: Option<&BTreeSet<String>>,
) -> bool {
    if let Some(range) = &options.range_filter {
        if transcript.seqid != range.seqid {
            return false;
        }
        if let Some(strand) = range.strand {
            if transcript.strand != strand {
                return false;
            }
        }
        if range.fully_within {
            if transcript.start < range.start || transcript.end > range.end {
                return false;
            }
        } else if transcript.start > range.end || transcript.end < range.start {
            return false;
        }
    }

    if let Some(id_filter) = &options.id_filter {
        let listed = id_set
            .expect("id set must exist when id filter is enabled")
            .contains(&transcript.id);
        match id_filter.mode {
            IdFilterMode::Include if !listed => return false,
            IdFilterMode::Exclude if listed => return false,
            IdFilterMode::Include | IdFilterMode::Exclude => {}
        }
    }

    if let Some(min_length) = options.min_length {
        if transcript.covlen() < min_length {
            return false;
        }
    }

    if let Some(max_intron) = options.max_intron {
        if transcript.max_intron_len() > max_intron {
            return false;
        }
    }

    if options.multi_exon_only && transcript.exons.len() <= 1 {
        return false;
    }

    if options.coding_only && !transcript.has_cds() {
        return false;
    }

    if options.noncoding_only && transcript.has_cds() {
        return false;
    }

    if options.no_pseudo && transcript.is_pseudo() {
        return false;
    }

    true
}

fn load_ids(path: &Path) -> Result<BTreeSet<String>, CompatError> {
    let text = fs::read_to_string(path).map_err(|_| {
        CompatError::new(
            format!("Error: cannot open input file {}!\n", path.display()),
            1,
        )
    })?;

    Ok(text
        .lines()
        .map(str::trim)
        .filter(|line| !line.is_empty())
        .map(str::to_owned)
        .collect())
}
