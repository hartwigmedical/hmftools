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
}

enum class VJGeneType
{
    IGHV,
    IGHJ,
    IGKV,
    IGKJ,
    IGLV,
    IGLJ,
    TRAV,
    TRAJ,
    TRBV,
    TRBJ,
    TRDV,
    TRDJ,
    TRGV,
    TRGJ;

    val isV: Boolean get() { return name[3] == 'V'; }
    val isJ: Boolean get() { return name[3] == 'J'; }
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
