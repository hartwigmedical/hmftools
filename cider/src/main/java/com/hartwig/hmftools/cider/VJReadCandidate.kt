package com.hartwig.hmftools.cider

import com.hartwig.hmftools.cider.genes.VJAnchorTemplate
import com.hartwig.hmftools.cider.genes.VJGeneType
import com.hartwig.hmftools.common.codon.Codons
import htsjdk.samtools.SAMRecord
import htsjdk.samtools.util.SequenceUtil

// We match reads to the genes, but they might not
// match perfectly. We just match by anchor
data class VJReadCandidate(
    val read: SAMRecord,
    val vjAnchorTemplates: List<VJAnchorTemplate>,
    val vjGeneType: VJGeneType,
    val templateAnchorSequence: String,
    val matchMethod: MatchMethod,
    val useReverseComplement: Boolean,
    val anchorOffsetStart: Int, // this is after reverse complement if needed, can be negative
    val anchorOffsetEnd: Int, // this is after reverse complement if needed, can be after sequence end
    val leftSoftClip: Int,
    val rightSoftClip: Int
)
{
    enum class MatchMethod
    {
        ALIGN, BLOSUM
    }

    var similarityScore: Int = Int.MIN_VALUE

    val readLength: Int get()
    {
        return read.readLength
    }

    // read sequence in the order that the gene is transcribed
    val readSequence: String get()
    {
        // now read out the sequence
        var readSeq: String = read.readString
        if (useReverseComplement)
            readSeq = SequenceUtil.reverseComplement(readSeq)
        return readSeq
    }

    val anchorSequence: String get()
    {
        if (anchorOffsetStart >= readLength || anchorOffsetEnd <= 0)
            return String()
        return readSequence.substring(Math.max(anchorOffsetStart, 0), Math.min(anchorOffsetEnd, readLength))
    }

    val anchorAA: String get()
    {
        return Codons.aminoAcidFromBases(anchorSequence)
    }
}
