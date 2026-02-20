package com.hartwig.hmftools.cider.genes

import com.hartwig.hmftools.common.genome.region.Strand
import com.hartwig.hmftools.common.region.ChrBaseRegion

// Represents a location inside the genome. Can be inside non primary assembly.
// posStart is 1 based and end is inclusive, consistent with reference genome convention
data class GenomicLocation(val chromosome: String,
                           val posStart: Int,
                           val posEnd: Int,
                           val strand: Strand,
                           val altAssemblyName: String? = null)
{
    init
    {
        require(posStart <= posEnd)
        require(altAssemblyName == null || altAssemblyName.isNotEmpty())
    }

    val inPrimaryAssembly: Boolean get() { return altAssemblyName == null }

    fun baseLength(): Int
    {
        return posEnd - posStart + 1
    }

    override fun toString(): String
    {
        return "${if (inPrimaryAssembly) "" else ("$altAssemblyName ") }${chromosome}:${posStart}-${posEnd}(${strand.asChar()})"
    }

    companion object
    {
        fun fromChrBaseRegionStrand(chrBaseRegion: ChrBaseRegion?, strand: Strand?): GenomicLocation?
        {
            return if (chrBaseRegion != null && strand != null)
                GenomicLocation(chrBaseRegion.chromosome(),
                    chrBaseRegion.start(),
                    chrBaseRegion.end(),
                    strand)
            else null
        }

        fun toChrBaseRegion(genomicLocation: GenomicLocation?): ChrBaseRegion?
        {
            return if (genomicLocation != null)
                ChrBaseRegion(genomicLocation.chromosome, genomicLocation.posStart, genomicLocation.posEnd)
            else null
        }
    }
}
