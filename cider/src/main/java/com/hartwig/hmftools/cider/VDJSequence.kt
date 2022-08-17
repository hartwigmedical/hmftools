package com.hartwig.hmftools.cider

import com.hartwig.hmftools.cider.layout.ReadLayout
import com.hartwig.hmftools.common.codon.Codons

interface VJAnchor
{
    enum class Type
    {
        V, J
    }

    val type: Type
    val geneType: VJGeneType
    val anchorBoundary: Int // end of V or start of J
    val matchMethod: String
}

data class VJAnchorByBlosum(
    override val type: VJAnchor.Type,
    override val anchorBoundary: Int,
    val templateAnchorSeq: String,
    val templateGenes: Collection<VJGene>,
    val similarityScore: Int
) : VJAnchor
{
    override val geneType: VJGeneType
    override val matchMethod: String get() = "blosum"

    init
    {
        val geneTypes: List<VJGeneType> = templateGenes.map({ o -> o.type }).distinct()
        if (geneTypes.size != 1)
            throw IllegalStateException("VJAnchorByBlosum: gene types(${geneTypes}) size != 1")
        geneType = geneTypes.first()
    }
}

data class VJAnchorByReadMatch(
    override val type: VJAnchor.Type,
    override val geneType: VJGeneType,
    override val anchorBoundary: Int,
    override val matchMethod: String
) : VJAnchor
{
}

class VDJSequence(
    val id: String,
    val layout: ReadLayout,
    val layoutSliceStart: Int,
    val layoutSliceEnd: Int,
    val vAnchor: VJAnchor,
    val jAnchor: VJAnchor)
{
    val numReads: Int get() = layout.reads.size

    init
    {
        //if (sequence.length != support.size)
          //  throw RuntimeException("VDJSequence: sequence.length != support.size")
    }

    val length: Int get()
    {
        return layoutSliceEnd - layoutSliceStart
    }

    val sequence: String get()
    {
        return layout.consensusSequence().substring(layoutSliceStart, layoutSliceEnd)
    }

    val sequenceFormatted: String get()
    {
        return CiderUtils.insertDashes(sequence, vAnchor.anchorBoundary, jAnchor.anchorBoundary)
    }

    val aminoAcidSequence: String get()
    {
        val codonAlignedSeq = sequence.drop(vAnchor.anchorBoundary % 3)
        return Codons.aminoAcidFromBases(codonAlignedSeq)
    }

    val aminoAcidSequenceFormatted: String get()
    {
        // we print the pre J anchor part first then the J anchor. This ensures that if J anchor is not aligned to codon
        // it will still get printed the way we want it
        val codonAlignedSeqBeforeJ = sequence.substring(0, jAnchor.anchorBoundary).drop(vAnchor.anchorBoundary % 3)

        return CiderUtils.insertDashes(Codons.aminoAcidFromBases(codonAlignedSeqBeforeJ), vAnchor.anchorBoundary / 3) + '-' +
                Codons.aminoAcidFromBases(jAnchorSequence)
    }

    val vAnchorLength: Int get()
    {
        return vAnchor.anchorBoundary
    }

    val jAnchorLength: Int get()
    {
        return length - jAnchor.anchorBoundary
    }

    val vAnchorSequence: String get()
    {
        return sequence.take(vAnchor.anchorBoundary)
    }

    val jAnchorSequence: String get()
    {
        return sequence.substring(jAnchor.anchorBoundary)
    }

    // does not include the C and W
    val cdr3SequenceShort: String get()
    {
        return sequence.substring(vAnchor.anchorBoundary, jAnchor.anchorBoundary)
    }

    val cdr3Sequence: String get()
    {
        return sequence.substring(Math.max(vAnchor.anchorBoundary - 3, 0), Math.min(jAnchor.anchorBoundary + 3, length))
    }

    val supportCounts: IntArray get()
    {
        return layout.highQualSupportCounts().sliceArray(layoutSliceStart until layoutSliceEnd)
    }

    val supportString: String get()
    {
        return CiderUtils.countsToString(supportCounts)
    }

    val supportMin: Int get() = supportCounts.minOrNull() ?: 0

    val cdr3SupportMin: Int get()
    {
        return supportCounts.slice(Math.max(vAnchor.anchorBoundary - 3, 0) until
                Math.min(jAnchor.anchorBoundary + 3, length)).minOrNull() ?: 0
    }

    val isInFrame: Boolean get()
    {
        return (jAnchor.anchorBoundary % 3) == 0
    }

    fun getSupportAt(index: Int) : Map.Entry<Char, Int>
    {
        return layout.getHighQualSequenceSupportAt(layoutSliceStart + index)
    }

}
