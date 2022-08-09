package com.hartwig.hmftools.cdr3

import com.hartwig.hmftools.common.codon.Codons
import com.hartwig.hmftools.common.genome.region.Strand
import com.hartwig.hmftools.common.utils.sv.ChrBaseRegion

data class GeneLocation(val chromosome: String, val posStart: Int, val posEnd: Int, val strand: Strand)
    : ChrBaseRegion(chromosome, posStart, posEnd)
{
    operator fun compareTo(other: GeneLocation): Int
    {
        val baseRegionCompare = super.compareTo(other)
        return if (baseRegionCompare == 0) strand.compareTo(other.strand) else baseRegionCompare
    }

    override fun toString(): String
    {
        return "${chromosome}:${posStart}-${posEnd}(${strand.asChar()})"
    }
}

enum class VJ
{
    V, J
}

enum class VJGeneType(val vj: VJ)
{
    IGHV(VJ.V),
    IGHJ(VJ.J),
    IGKV(VJ.V),
    IGKJ(VJ.J),
    IGLV(VJ.V),
    IGLJ(VJ.J),
    TRAV(VJ.V),
    TRAJ(VJ.J),
    TRBV(VJ.V),
    TRBJ(VJ.J),
    TRDV(VJ.V),
    TRDJ(VJ.J),
    TRGV(VJ.V),
    TRGJ(VJ.J);
}

// store the anchor location also whether it is V or J
data class VJAnchorReferenceLocation(val vj: VJ, val geneLocation: GeneLocation)
{
    val chromosome: String get() = geneLocation.chromosome
    val start: Int get() = geneLocation.posStart
    val end: Int get() = geneLocation.posEnd
    val strand: Strand get() = geneLocation.strand

    fun baseLength() : Int
    {
        return geneLocation.baseLength()
    }

    // get the reference position of the end of the anchor
    // this is the end of the C codon for V and the start of the W / F codon for J
    fun anchorEndReferencePosition() : Int
    {
        if (vj == VJ.V && geneLocation.strand == Strand.FORWARD ||
            vj == VJ.J && geneLocation.strand == Strand.REVERSE)
        {
            return geneLocation.posEnd
        }
        else
        {
            return geneLocation.posStart
        }
    }
}

data class VJGene
    (
    val id: String,
    val name: String,
    val allele: String,
    val geneLocation: GeneLocation?,
    val sequence: String,
    val anchorSequence: String,
    val anchorLocation: GeneLocation?
)
{
    val type: VJGeneType = VJGeneType.valueOf(name.take(4))
    val anchorAminoAcidSequence: String = Codons.aminoAcidFromBases(anchorSequence)
    val chromosome: String? get() { return geneLocation?.chromosome() }
    val startPosition: Int get() { return geneLocation?.start() ?: -1 }
    val endPosition: Int get() { return geneLocation?.end() ?: -1 }
    val strand: Strand? get() { return geneLocation?.strand }
}

fun getConservedAA(vjGeneType: VJGeneType) : Char
{
    return when (vjGeneType)
    {
        VJGeneType.IGHV -> 'C'
        VJGeneType.IGHJ -> 'W'
        VJGeneType.IGKV -> 'C'
        VJGeneType.IGKJ -> 'F'
        VJGeneType.IGLV -> 'C'
        VJGeneType.IGLJ -> 'F'
        VJGeneType.TRAV -> 'C'
        VJGeneType.TRAJ -> 'F'
        VJGeneType.TRBV -> 'C'
        VJGeneType.TRBJ -> 'F'
        VJGeneType.TRDV -> 'C'
        VJGeneType.TRDJ -> 'F'
        VJGeneType.TRGV -> 'C'
        VJGeneType.TRGJ -> 'F'
    }
}
