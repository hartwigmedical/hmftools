package com.hartwig.hmftools.cider

import com.hartwig.hmftools.common.codon.Codons
import com.hartwig.hmftools.common.genome.region.Strand
import com.hartwig.hmftools.common.utils.sv.ChrBaseRegion

data class GenomeRegionStrand(val chromosome: String, val posStart: Int, val posEnd: Int, val strand: Strand)
    : ChrBaseRegion(chromosome, posStart, posEnd)
{
    operator fun compareTo(other: GenomeRegionStrand): Int
    {
        val baseRegionCompare = super.compareTo(other)
        return if (baseRegionCompare == 0) strand.compareTo(other.strand) else baseRegionCompare
    }

    override fun toString(): String
    {
        return "${chromosome}:${posStart}-${posEnd}(${strand.asChar()})"
    }
}

data class VJAnchorTemplate
    (
    val id: String,
    val name: String,
    val allele: String,
    val geneLocation: GenomeRegionStrand?,
    val sequence: String,
    val anchorSequence: String,
    val anchorLocation: GenomeRegionStrand?
)
{
    val type: VJGeneType = VJGeneType.valueOf(name.take(4))

    val vj: VJ get() { return type.vj }
    val anchorAminoAcidSequence: String = Codons.aminoAcidFromBases(anchorSequence)
    val chromosome: String? get() { return geneLocation?.chromosome() }
    val startPosition: Int get() { return geneLocation?.start() ?: -1 }
    val endPosition: Int get() { return geneLocation?.end() ?: -1 }
    val strand: Strand? get() { return geneLocation?.strand }
}

// store the anchor location and also the type of the gene segment
data class VJAnchorGenomeLocation(val vjGeneType: VJGeneType, val genomeLocation: GenomeRegionStrand)
{
    val vj: VJ get() = vjGeneType.vj
    val chromosome: String get() = genomeLocation.chromosome
    val start: Int get() = genomeLocation.posStart
    val end: Int get() = genomeLocation.posEnd
    val strand: Strand get() = genomeLocation.strand

    fun baseLength() : Int
    {
        return genomeLocation.baseLength()
    }

    // get the reference position of the end of the anchor
    // this is the end of the C codon for V and the start of the W / F codon for J
    fun anchorBoundaryReferencePosition() : Int
    {
        return if (vjGeneType.vj == VJ.V && genomeLocation.strand == Strand.FORWARD ||
            vjGeneType.vj == VJ.J && genomeLocation.strand == Strand.REVERSE)
            {
                genomeLocation.posEnd
            }
            else
            {
                genomeLocation.posStart
            }
    }
}
