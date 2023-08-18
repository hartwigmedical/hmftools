package com.hartwig.hmftools.cider.genes

import com.hartwig.hmftools.common.genome.region.Strand

// Represents a location inside the genome. Can be inside non primary assembly.
// posStart is 1 based and end is inclusive, consistent with reference genome convention
data class GenomicLocation(val isPrimaryAssembly: Boolean,
                           val assemblyName: String,
                           val chromosome: String,
                           val posStart: Int,
                           val posEnd: Int,
                           val strand: Strand)
{
    fun baseLength(): Int
    {
        return posEnd - posStart + 1
    }

    /*
    operator fun compareTo(other: GenomicLocation): Int
    {
        val baseRegionCompare = super.compareTo(other)
        return if (baseRegionCompare == 0) strand.compareTo(other.strand) else baseRegionCompare
    }*/

    override fun toString(): String
    {
        return "${if (isPrimaryAssembly) "" else "alt "}${chromosome}:${posStart}-${posEnd}(${strand.asChar()})"
    }
}
