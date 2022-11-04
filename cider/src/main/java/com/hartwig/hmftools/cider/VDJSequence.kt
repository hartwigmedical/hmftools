package com.hartwig.hmftools.cider

import com.hartwig.hmftools.cider.layout.ReadLayout
import com.hartwig.hmftools.common.codon.Codons
import org.eclipse.collections.api.collection.ImmutableCollection

interface VJAnchor
{
    val vj: VJ
    val geneType: VJGeneType
    val anchorBoundary: Int // last base of V + 1 or first base of J
    val matchMethod: String
    val templateAnchorSeq: String
}

data class VJAnchorByBlosum(
    override val vj: VJ,
    override val geneType: VJGeneType,
    override val anchorBoundary: Int,
    override val templateAnchorSeq: String,
    val templateGenes: ImmutableCollection<VJAnchorTemplate>,
    val similarityScore: Int
) : VJAnchor
{
    override val matchMethod: String get() = "blosum"

    init
    {
        val geneTypeList: List<VJGeneType> = templateGenes.map({ o -> o.type }).distinct()
        if (geneTypeList.size != 1)
            throw IllegalStateException("VJAnchorByBlosum: gene types(${geneTypeList}) size != 1")
        if (geneTypeList[0] != geneType)
            throw IllegalStateException("VJAnchorByBlosum: gene type mismatch: ${geneType} != ${geneTypeList[0]}")
    }
}

data class VJAnchorByReadMatch(
    override val vj: VJ,
    override val geneType: VJGeneType,
    override val anchorBoundary: Int,
    override val matchMethod: String,
    override val templateAnchorSeq: String,
    val numReads: Int
) : VJAnchor
{
}

class VDJSequence(
    val layout: ReadLayout,
    val layoutSliceStart: Int,
    val layoutSliceEnd: Int,
    val vAnchor: VJAnchor?,
    val jAnchor: VJAnchor?)
{
    val numReads: Int get() = layout.reads.size

    val length: Int get()
    {
        return layoutSliceEnd - layoutSliceStart
    }

    val sequence: String get()
    {
        return layout.consensusSequence().substring(layoutSliceStart, layoutSliceEnd)
    }

    val isFullyRearranged: Boolean get()
    {
        return vAnchor != null && jAnchor != null
    }

    val vAnchorBoundary: Int? get()
    {
        return vAnchor?.anchorBoundary
    }

    val jAnchorBoundary: Int? get()
    {
        return jAnchor?.anchorBoundary
    }

    val sequenceFormatted: String get()
    {
        return CiderUtils.insertDashes(sequence, vAnchorBoundary ?: 0, jAnchorBoundary ?: length)
    }

    val aminoAcidSequence: String get()
    {
        val codonAlignedSeq = sequence.drop((vAnchorBoundary ?: jAnchorBoundary ?: 0) % 3)
        return Codons.aminoAcidFromBases(codonAlignedSeq)
    }

    val aminoAcidSequenceFormatted: String get()
    {
        if (vAnchor == null)
        {
            if (jAnchor == null)
                return String()

            // if there is no v anchor, we have to align to j anchor
            // also align to codon boundary for the J anchor
            val codonAlignedSeq = sequence.substring(jAnchor.anchorBoundary % 3)
            val dashPos = jAnchor.anchorBoundary / 3
            return CiderUtils.insertDashes(Codons.aminoAcidFromBases(codonAlignedSeq), dashPos)
        }
        else if (jAnchor == null)
        {
            // if no j anchor, we just need to align to V anchor
            val codonAlignedSeq = sequence.substring(vAnchor.anchorBoundary % 3)
            val dashPos = vAnchor.anchorBoundary / 3
            return CiderUtils.insertDashes(Codons.aminoAcidFromBases(codonAlignedSeq), dashPos)
        }

        // we print the pre J anchor part first then the J anchor. This ensures that if J anchor is not aligned to codon
        // it will still get printed the way we want it
        val codonAlignedSeqBeforeJ = sequence.substring(vAnchor.anchorBoundary % 3, jAnchor.anchorBoundary)
        val dashPos = vAnchor.anchorBoundary / 3

        return CiderUtils.insertDashes(Codons.aminoAcidFromBases(codonAlignedSeqBeforeJ), dashPos) +
                '-' +
                Codons.aminoAcidFromBases(jAnchorSequence)
    }

    val vAnchorLength: Int get()
    {
        return vAnchor?.anchorBoundary ?: 0
    }

    val jAnchorLength: Int get()
    {
        return length - (jAnchor?.anchorBoundary ?: 0)
    }

    val vAnchorSequence: String get()
    {
        return sequence.take(vAnchorBoundary ?: 0)
    }

    val jAnchorSequence: String get()
    {
        return if (jAnchor == null) String() else sequence.substring(jAnchorBoundary!!)
    }

    // does not include the C and W
    val cdr3SequenceShort: String get()
    {
        if (vAnchorBoundary != null)
        {
            return sequence.substring(vAnchorBoundary ?: 0, jAnchorBoundary ?: length)
        }
        if (jAnchorBoundary != null)
        {
            // we need to aligned it to the jAnchor boundary, avoid it being out of frame
            return sequence.substring(jAnchorBoundary!! % 3, jAnchorBoundary!!)
        }
        // shouldn't get here
        return sequence
    }

    private val cdr3Start: Int get()
    {
        if (vAnchorBoundary != null)
        {
            return Math.max(vAnchorBoundary!! - 3, 0)
        }
        if (jAnchorBoundary != null)
        {
            // we need to aligned it to the jAnchor boundary, avoid it being out of frame
            return jAnchorBoundary!! % 3
        }
        // shouldn't get here
        return 0
    }

    private val cdr3End: Int get()
    {
        return Math.min((jAnchorBoundary ?: return length) + 3, length)
    }

    val cdr3Sequence: String get()
    {
        return sequence.substring(cdr3Start, cdr3End)
    }

    val supportCounts: IntArray get()
    {
        return layout.highQualSupportCounts().sliceArray(layoutSliceStart until layoutSliceEnd)
    }

    val supportString: String get()
    {
        return CiderUtils.countsToString(supportCounts)
    }

    val cdr3SupportMin: Int get()
    {
        return supportCounts.slice(cdr3Start until cdr3End).minOrNull() ?: 0
    }

    val isInFrame: Boolean get()
    {
        return ((jAnchorBoundary ?: 0) % 3) == 0
    }

    fun getSupportAt(index: Int) : Map.Entry<Char, Int>
    {
        return layout.getHighQualSequenceSupportAt(layoutSliceStart + index)
    }

}
