use std::cmp::Ordering;
use std::collections::BTreeMap;
use std::fs;

use crate::compat::CompatError;
use crate::model::{Annotation, Gene, Transcript};
use crate::options::RefSortOrder;

pub fn sort_transcripts_for_fasta(
    transcripts: &[Transcript],
    order: &RefSortOrder,
    ref_order: &[String],
) -> Result<Vec<Transcript>, CompatError> {
    let ranked_refs = rank_refs(order, ref_order)?;
    let mut result = transcripts.to_vec();
    gclib_qsort(&mut result, |left, right| {
        compare_transcript_by_loc(left, right, &ranked_refs)
    });
    Ok(result)
}

pub fn sort_annotation(
    annotation: &mut Annotation,
    order: &RefSortOrder,
    keep_genes: bool,
) -> Result<(), CompatError> {
    let ref_order = ordered_refs(order, &annotation.ref_order)?;

    let ranks = ref_order
        .iter()
        .enumerate()
        .map(|(index, seqid)| (seqid.clone(), index))
        .collect::<BTreeMap<_, _>>();
    let mut records = sort_records_in_original_order(annotation, keep_genes);
    gclib_qsort(&mut records, |left, right| {
        compare_record_by_loc(left, right, annotation, &ranks)
    });

    let mut genes = Vec::with_capacity(annotation.genes.len());
    let mut transcripts = Vec::with_capacity(annotation.transcripts.len());
    for record in records {
        match record {
            SortRecord::Gene(index) => gclib_sorted_insert(
                &mut genes,
                annotation.genes[index].clone(),
                |left, right| compare_gene_by_loc(left, right, &ranks),
            ),
            SortRecord::Transcript(index) => gclib_sorted_insert(
                &mut transcripts,
                annotation.transcripts[index].clone(),
                |left, right| compare_transcript_by_loc(left, right, &ranks),
            ),
        }
    }
    if keep_genes {
        annotation.genes = genes;
    } else {
        gclib_qsort(&mut annotation.genes, |left, right| {
            compare_gene_by_loc(left, right, &ranks)
        });
    }
    annotation.transcripts = transcripts;
    annotation.ref_order = ref_order;

    Ok(())
}

fn ordered_refs(order: &RefSortOrder, ref_order: &[String]) -> Result<Vec<String>, CompatError> {
    Ok(match order {
        RefSortOrder::Input => ref_order.to_vec(),
        RefSortOrder::Alpha => {
            let mut refs = ref_order.to_vec();
            refs.sort();
            refs
        }
        RefSortOrder::List(path) => {
            let mut refs = parse_ref_list(&fs::read_to_string(path).map_err(|_| {
                CompatError::new(
                    format!("Error: could not open file {} !\n", path.display()),
                    1,
                )
            })?);
            for seqid in ref_order {
                if !refs.iter().any(|candidate| candidate == seqid) {
                    refs.push(seqid.clone());
                }
            }
            refs
        }
    })
}

fn rank_refs(
    order: &RefSortOrder,
    ref_order: &[String],
) -> Result<BTreeMap<String, usize>, CompatError> {
    Ok(ordered_refs(order, ref_order)?
        .iter()
        .enumerate()
        .map(|(index, seqid)| (seqid.clone(), index))
        .collect::<BTreeMap<_, _>>())
}

#[derive(Clone, Copy)]
enum SortRecord {
    Gene(usize),
    Transcript(usize),
}

fn sort_records_in_original_order(annotation: &Annotation, keep_genes: bool) -> Vec<SortRecord> {
    let mut records = Vec::with_capacity(annotation.genes.len() + annotation.transcripts.len());
    if keep_genes {
        records.extend((0..annotation.genes.len()).map(SortRecord::Gene));
    }
    records.extend((0..annotation.transcripts.len()).map(SortRecord::Transcript));
    records.sort_by_key(|record| match record {
        SortRecord::Gene(index) => annotation.genes[*index].order,
        SortRecord::Transcript(index) => annotation.transcripts[*index].order,
    });
    records
}

fn compare_record_by_loc(
    left: &SortRecord,
    right: &SortRecord,
    annotation: &Annotation,
    ranks: &BTreeMap<String, usize>,
) -> Ordering {
    match (left, right) {
        (SortRecord::Gene(left), SortRecord::Gene(right)) => {
            compare_gene_by_loc(&annotation.genes[*left], &annotation.genes[*right], ranks)
        }
        (SortRecord::Transcript(left), SortRecord::Transcript(right)) => compare_transcript_by_loc(
            &annotation.transcripts[*left],
            &annotation.transcripts[*right],
            ranks,
        ),
        _ => record_key(left, annotation, ranks).cmp(&record_key(right, annotation, ranks)),
    }
}

fn record_key<'a>(
    record: &SortRecord,
    annotation: &'a Annotation,
    ranks: &BTreeMap<String, usize>,
) -> (usize, u64, u8, u64, &'a str) {
    match record {
        SortRecord::Gene(index) => {
            let gene = &annotation.genes[*index];
            (
                ranks
                    .get(gene.seqid.as_str())
                    .copied()
                    .unwrap_or(usize::MAX),
                gene.start,
                0,
                gene.end,
                gene.id.as_str(),
            )
        }
        SortRecord::Transcript(index) => {
            let transcript = &annotation.transcripts[*index];
            (
                ranks
                    .get(transcript.seqid.as_str())
                    .copied()
                    .unwrap_or(usize::MAX),
                transcript.start,
                feature_level(transcript.gene_id.is_some()),
                transcript.end,
                transcript.id.as_str(),
            )
        }
    }
}

fn compare_gene_by_loc(left: &Gene, right: &Gene, ranks: &BTreeMap<String, usize>) -> Ordering {
    (
        ranks
            .get(left.seqid.as_str())
            .copied()
            .unwrap_or(usize::MAX),
        left.start,
        0,
        left.end,
        &left.id,
    )
        .cmp(&(
            ranks
                .get(right.seqid.as_str())
                .copied()
                .unwrap_or(usize::MAX),
            right.start,
            0,
            right.end,
            &right.id,
        ))
}

fn compare_transcript_by_loc(
    left: &Transcript,
    right: &Transcript,
    ranks: &BTreeMap<String, usize>,
) -> Ordering {
    (
        ranks
            .get(left.seqid.as_str())
            .copied()
            .unwrap_or(usize::MAX),
        left.start,
        feature_level(left.gene_id.is_some()),
        left.end,
        &left.id,
    )
        .cmp(&(
            ranks
                .get(right.seqid.as_str())
                .copied()
                .unwrap_or(usize::MAX),
            right.start,
            feature_level(right.gene_id.is_some()),
            right.end,
            &right.id,
        ))
}

fn feature_level(has_parent: bool) -> u8 {
    if has_parent {
        1
    } else {
        0
    }
}

fn gclib_qsort<T, F>(items: &mut [T], compare: F)
where
    T: Clone,
    F: Fn(&T, &T) -> Ordering,
{
    if items.len() <= 1 {
        return;
    }
    gclib_qsort_range(items, 0, items.len() - 1, &compare);
}

fn gclib_qsort_range<T, F>(items: &mut [T], mut left: usize, right: usize, compare: &F)
where
    T: Clone,
    F: Fn(&T, &T) -> Ordering,
{
    loop {
        let mut i = left;
        let mut j = right;
        let pivot = items[left + ((right - left) >> 1)].clone();
        loop {
            while compare(&items[i], &pivot).is_lt() {
                i += 1;
            }
            while compare(&items[j], &pivot).is_gt() {
                j -= 1;
            }
            if i <= j {
                items.swap(i, j);
                i += 1;
                if j == 0 {
                    break;
                }
                j -= 1;
            }
            if i > j {
                break;
            }
        }
        if left < j {
            gclib_qsort_range(items, left, j, compare);
        }
        left = i;
        if i >= right {
            break;
        }
    }
}

fn gclib_sorted_insert<T, F>(items: &mut Vec<T>, item: T, compare: F)
where
    F: Fn(&T, &T) -> Ordering,
{
    let index = gclib_insert_index(items, &item, &compare);
    items.insert(index, item);
}

fn gclib_insert_index<T, F>(items: &[T], item: &T, compare: &F) -> usize
where
    F: Fn(&T, &T) -> Ordering,
{
    if items.is_empty() {
        return 0;
    }
    if compare(&items[0], item).is_gt() {
        return 0;
    }
    if compare(item, items.last().expect("checked non-empty")).is_gt() {
        return items.len();
    }

    let mut low = 0usize;
    let mut high = items.len() - 1;
    while low <= high {
        let index = low + ((high - low) >> 1);
        match compare(&items[index], item) {
            Ordering::Less => low = index + 1,
            Ordering::Equal => return index,
            Ordering::Greater => {
                if index == 0 {
                    return 0;
                }
                high = index - 1;
            }
        }
    }
    low
}

fn parse_ref_list(text: &str) -> Vec<String> {
    text.split(|ch: char| ch.is_whitespace() || matches!(ch, ',' | ';'))
        .filter(|part| !part.is_empty())
        .map(str::to_owned)
        .collect()
}
