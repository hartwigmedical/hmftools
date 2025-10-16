package com.hartwig.hmftools.cider.genes

import com.hartwig.hmftools.common.cider.IgTcrGene

fun IgTcrGene.genomicLocation(): GenomicLocation?
{
    return GenomicLocation.fromNullableFields(contigName, genePosition, geneStrand)
}

fun IgTcrGene.anchorGenomicLocation(): GenomicLocation?
{
    return GenomicLocation.fromNullableFields(contigName, anchorPosition, geneStrand)
}
