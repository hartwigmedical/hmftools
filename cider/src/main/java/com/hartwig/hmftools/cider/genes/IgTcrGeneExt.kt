package com.hartwig.hmftools.cider.genes

import com.hartwig.hmftools.common.cider.IgTcrGene

fun IgTcrGene.contig() = if (contigName == null) null else Contig.fromName(contigName!!)

fun IgTcrGene.genomicLocation(): GenomicLocation?
{
    return GenomicLocation.fromNullableFields(contig(), genePosition, geneStrand)
}

fun IgTcrGene.anchorGenomicLocation(): GenomicLocation?
{
    return GenomicLocation.fromNullableFields(contig(), anchorPosition, geneStrand)
}
