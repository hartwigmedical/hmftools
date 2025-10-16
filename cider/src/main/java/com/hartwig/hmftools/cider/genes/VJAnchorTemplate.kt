package com.hartwig.hmftools.cider.genes

import com.hartwig.hmftools.cider.VJ
import com.hartwig.hmftools.cider.VJGeneType
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
    val strand: Strand? get() { return geneLocation?.strand }
}

// store the anchor location and also the type of the gene segment
data class VJAnchorGenomeLocation(val vjGeneType: VJGeneType, val genomeLocation: GenomicLocation)
{
    val vj: VJ get() = vjGeneType.vj
    val contig: Contig get() = genomeLocation.contig
    val start: Int get() = genomeLocation.position.start()
    val end: Int get() = genomeLocation.position.end()
    val strand: Strand get() = genomeLocation.strand

    fun baseLength() : Int
    {
        return genomeLocation.position.baseLength()
    }

    // get the reference position of the end of the anchor
    // this is the end of the C codon for V and the start of the W / F codon for J
    fun anchorBoundaryReferencePosition() : Int
    {
        return if (vjGeneType.vj == VJ.V && genomeLocation.strand == Strand.FORWARD ||
            vjGeneType.vj == VJ.J && genomeLocation.strand == Strand.REVERSE)
            {
                genomeLocation.position.end()
            }
            else
            {
                genomeLocation.position.start()
            }
    }
}
