package com.hartwig.hmftools.cdr3

import com.hartwig.hmftools.cdr3.layout.ReadLayout
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
    val layoutStart: Int,
    val layoutEnd: Int,
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
        return layoutEnd - layoutStart
    }

    val sequence: String get()
    {
        return layout.consensusSequence().substring(layoutStart, layoutEnd)
    }

    val sequenceFormatted: String get()
    {
        return Cdr3Utils.insertDashes(sequence, vAnchor.anchorBoundary, jAnchor.anchorBoundary)
    }

    val aminoAcidSequence: String get()
    {
        return Codons.aminoAcidFromBases(sequence)
    }

    val aminoAcidSequenceFormatted: String get()
    {
        return Cdr3Utils.insertDashes(aminoAcidSequence, vAnchor.anchorBoundary / 3, jAnchor.anchorBoundary / 3)
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
    val dSequenceShort: String get()
    {
        return sequence.substring(vAnchor.anchorBoundary, jAnchor.anchorBoundary)
    }

    val dSequence: String get()
    {
        return sequence.substring(vAnchor.anchorBoundary - 3, jAnchor.anchorBoundary + 3)
    }

    val supportCounts: IntArray get()
    {
        return layout.highQualSupportCounts().sliceArray(layoutStart until layoutEnd)
    }

    val supportString: String get()
    {
        return Cdr3Utils.countsToString(supportCounts)
    }

    val supportMin: Int get() = supportCounts.minOrNull() ?: 0

    val cdr3SupportMin: Int get()
    {
        return supportCounts.slice(vAnchor.anchorBoundary - 3 until jAnchor.anchorBoundary + 3).minOrNull() ?: 0
    }

    fun getSupportAt(index: Int) : Map.Entry<Char, Int>
    {
        return layout.getHighQualSequenceSupportAt(layoutStart + index)
    }

}
