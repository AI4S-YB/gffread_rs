use std::io::{self, Write};

use crate::model::Locus;

pub fn write_locus<W: Write>(out: &mut W, locus: &Locus) -> io::Result<()> {
    writeln!(
        out,
        "{}\tgffcl\tlocus\t{}\t{}\t.\t{}\t.\tID={};genes={};geneIDs={};transcripts={}",
        locus.seqid,
        locus.start,
        locus.end,
        locus.strand,
        locus.id,
        locus.genes.join(","),
        locus.gene_ids.join(","),
        locus.transcripts.join(",")
    )
}
