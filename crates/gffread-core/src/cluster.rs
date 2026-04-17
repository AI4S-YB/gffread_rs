use std::collections::BTreeSet;

use crate::model::{Annotation, Locus, Transcript};
use crate::options::ClusterOptions;

pub fn apply_clustering(
    annotation: &Annotation,
    options: &ClusterOptions,
) -> (Annotation, Vec<Locus>) {
    if !options.merge {
        return (annotation.clone(), Vec::new());
    }

    let mut transcripts = annotation.transcripts.clone();
    transcripts.sort_by(|left, right| {
        (&left.seqid, left.start, left.end, &left.id).cmp(&(
            &right.seqid,
            right.start,
            right.end,
            &right.id,
        ))
    });

    let mut loci = Vec::new();
    let mut result = Vec::new();
    let mut index = 0usize;
    let mut locus_index = 1usize;

    while index < transcripts.len() {
        let mut group = vec![transcripts[index].clone()];
        let mut end = transcripts[index].end;
        let seqid = transcripts[index].seqid.clone();

        let mut cursor = index + 1;
        while cursor < transcripts.len()
            && transcripts[cursor].seqid == seqid
            && transcripts[cursor].start <= end
        {
            end = end.max(transcripts[cursor].end);
            group.push(transcripts[cursor].clone());
            cursor += 1;
        }

        let locus_id = format!("RLOC_{locus_index:08}");
        locus_index += 1;
        let kept = if options.cluster_only {
            group.clone()
        } else {
            collapse_redundancy(group.clone(), options)
        };

        let locus = make_locus(&group, &kept, &locus_id);

        for mut transcript in kept {
            transcript.locus = Some(locus_id.clone());
            result.push(transcript);
        }

        loci.push(locus);
        index = cursor;
    }

    (
        Annotation {
            transcripts: result,
        },
        loci,
    )
}

fn make_locus(group: &[Transcript], kept: &[Transcript], locus_id: &str) -> Locus {
    let mut genes = BTreeSet::new();
    let mut gene_ids = BTreeSet::new();
    for transcript in group {
        if let Some(gene) = &transcript.gene_name {
            genes.insert(gene.clone());
        }
        if let Some(gene_id) = &transcript.gene_id {
            gene_ids.insert(gene_id.clone());
        }
    }

    let strand = if group.iter().all(|transcript| transcript.strand == '+') {
        '+'
    } else if group.iter().all(|transcript| transcript.strand == '-') {
        '-'
    } else {
        '.'
    };

    Locus {
        seqid: group[0].seqid.clone(),
        start: group
            .iter()
            .map(|transcript| transcript.start)
            .min()
            .unwrap_or(0),
        end: group
            .iter()
            .map(|transcript| transcript.end)
            .max()
            .unwrap_or(0),
        strand,
        id: locus_id.to_owned(),
        genes: genes.into_iter().collect(),
        gene_ids: gene_ids.into_iter().collect(),
        transcripts: kept
            .iter()
            .map(|transcript| transcript.id.clone())
            .collect(),
    }
}

fn collapse_redundancy(group: Vec<Transcript>, options: &ClusterOptions) -> Vec<Transcript> {
    let mut kept: Vec<Transcript> = Vec::new();

    for transcript in group {
        if let Some(existing) = kept
            .iter_mut()
            .find(|existing| redundant(existing, &transcript, options))
        {
            if transcript.covlen() > existing.covlen() {
                *existing = transcript;
            }
        } else {
            kept.push(transcript);
        }
    }

    kept
}

fn redundant(left: &Transcript, right: &Transcript, options: &ClusterOptions) -> bool {
    if exact_intron_chain(left, right) {
        return options.relax_boundary_containment
            || contains_span(left, right)
            || contains_span(right, left);
    }

    options.collapse_contained && contained_intron_chain(left, right)
}

fn contains_span(container: &Transcript, contained: &Transcript) -> bool {
    container.start <= contained.start && container.end >= contained.end
}

fn exact_intron_chain(left: &Transcript, right: &Transcript) -> bool {
    if left.exons.len() != right.exons.len() {
        return false;
    }
    if left.exons.len() <= 1 {
        return contains_span(left, right) || contains_span(right, left);
    }

    left.exons
        .windows(2)
        .zip(right.exons.windows(2))
        .all(|(left_pair, right_pair)| {
            left_pair[0].end == right_pair[0].end && left_pair[1].start == right_pair[1].start
        })
}

fn contained_intron_chain(left: &Transcript, right: &Transcript) -> bool {
    if left.exons.len() == right.exons.len() {
        return exact_intron_chain(left, right);
    }

    let (shorter, longer) = if left.exons.len() < right.exons.len() {
        (left, right)
    } else {
        (right, left)
    };

    if !contains_span(longer, shorter) {
        return false;
    }

    let shorter_introns = introns(shorter);
    let longer_introns = introns(longer);
    shorter_introns
        .iter()
        .all(|intron| longer_introns.iter().any(|candidate| candidate == intron))
}

fn introns(transcript: &Transcript) -> Vec<(u64, u64)> {
    transcript
        .exons
        .windows(2)
        .map(|pair| (pair[0].end, pair[1].start))
        .collect()
}
