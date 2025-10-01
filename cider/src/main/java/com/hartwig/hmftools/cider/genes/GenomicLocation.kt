package com.hartwig.hmftools.cider.genes

import com.hartwig.hmftools.common.genome.region.Strand
import com.hartwig.hmftools.common.region.BaseRegion

data class GenomicLocation(
    val contig: Contig,
    val position: BaseRegion,
    val strand: Strand,
)
{
    init
    {
        require(position.start() <= position.end())
    }

    override fun toString(): String
    {
        return "${contig}:${position}(${strand.asChar()})"
    }

    companion object
    {
        fun fromNullableFields(contig: Contig?, position: BaseRegion?, strand: Strand?): GenomicLocation?
        {
            return if (contig != null && position != null && strand != null)
            {
                GenomicLocation(contig, position, strand)
            }
            else
            {
                null
            }
        }
    }
}
