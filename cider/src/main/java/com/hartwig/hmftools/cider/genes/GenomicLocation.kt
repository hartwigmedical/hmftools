package com.hartwig.hmftools.cider.genes

import com.hartwig.hmftools.common.cider.IgTcrGene
import com.hartwig.hmftools.common.genome.region.Strand
import com.hartwig.hmftools.common.region.ChrBaseRegion

// Represents a location inside the genome. Can be inside non primary assembly.
// posStart is 1 based and end is inclusive, consistent with reference genome convention
data class GenomicLocation(val bases: ChrBaseRegion,
                           val strand: Strand,
                           val inPrimaryAssembly: Boolean = true)
{
    init
    {
        require(bases.start() <= bases.end())
    }

    override fun toString(): String
    {
        val chr = if (inPrimaryAssembly) bases.chromosome() else "${bases.chromosome()}*"
        return "${chr}:${bases.start()}-${bases.end()}(${strand.asChar()})"
    }

    companion object
    {
        fun fromNullableFields(region: ChrBaseRegion?, strand: Strand?, inPrimaryAssembly: Boolean?): GenomicLocation?
        {
            return if (region != null && strand != null && inPrimaryAssembly != null)
            {
                GenomicLocation(region, strand, inPrimaryAssembly)
            }
            else
            {
                null
            }
        }
    }
}

fun IgTcrGene.genomicLocation(): GenomicLocation?
{
    return GenomicLocation.fromNullableFields(geneLocation, geneStrand, inPrimaryAssembly)
}

fun IgTcrGene.anchorGenomicLocation(): GenomicLocation?
{
    return GenomicLocation.fromNullableFields(anchorLocation, geneStrand, inPrimaryAssembly)
}
