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
    //val startPosition: Int get() { return geneLocation?.start() ?: -1 }
    //val endPosition: Int get() { return geneLocation?.end() ?: -1 }
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
