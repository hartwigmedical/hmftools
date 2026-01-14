package com.hartwig.hmftools.cider

import com.hartwig.hmftools.cider.genes.GenomicLocation
import com.hartwig.hmftools.common.codon.Codons
import com.hartwig.hmftools.common.genome.region.Strand

data class VJAnchorTemplate
    (
    val type: VJGeneType,
    val geneName: String, // IGHV1-45
    val allele: String, // 01
    val geneLocation: GenomicLocation?,
    val anchorSequence: String,
    val anchorLocation: GenomicLocation?
)
{
    init
    {
        require(anchorSequence.isNotEmpty())
        require(anchorLocation == null || anchorLocation.inPrimaryAssembly)
    }

    val vj: VJ get() { return type.vj }
    val anchorAminoAcidSequence: String = Codons.aminoAcidFromBases(anchorSequence)
    val chromosome: String? get() { return geneLocation?.chromosome }
    val strand: Strand? get() { return geneLocation?.strand }
}

// store the anchor location and also the type of the gene segment
data class VJAnchorGenomeLocation(val vjGeneType: VJGeneType, val genomeLocation: GenomicLocation)
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

    // get the side of the anchor towards the CDR3 sequence
    // this is the end of the C codon for V and the start of the W / F codon for J
    fun anchorBoundarySide() : Int
    {
        val rightSide = vjGeneType.vj == VJ.V && genomeLocation.strand == Strand.FORWARD ||
                vjGeneType.vj == VJ.J && genomeLocation.strand == Strand.REVERSE
        return if (rightSide) 1 else -1
    }

    fun anchorBoundaryReferencePosition() : Int
    {
        return if (anchorBoundarySide() == 1) end else start
    }
}
