package com.hartwig.hmftools.cdr3

import com.hartwig.hmftools.common.codon.Codons
import com.hartwig.hmftools.common.genome.region.GenomeRegion
import htsjdk.samtools.SAMRecord
import htsjdk.samtools.util.SequenceUtil

// this one includes supplementary etc
data class VJReadRecordKey(
    val readName: String,
    val firstOfPair: Boolean,
    val alignment: GenomeRegion?
)

// We match reads to the genes, but they might not
// match perfectly. We just match by anchor
data class VJReadCandidate(
    val read: SAMRecord,
    val vjGenes: List<VJGene>,
    val vjGeneType: VJGeneType,
    val anchorMatchType: AnchorMatchType,
    val useReverseComplement: Boolean,
    val anchorOffsetStart: Int, // this is after reverse complement if needed, can be negative
    val anchorOffsetEnd: Int, // this is after reverse complement if needed, can be after sequence end
    val anchorLocation: GeneLocation?, // could be null if not mapped
    val leftSoftClip: Int,
    val rightSoftClip: Int
)
{
    enum class AnchorMatchType
    {
        ALIGN, EXACT
    }

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

    // base qualities in the order that the gene is transcribed
    val baseQualities: ByteArray get()
    {
        var bq = read.baseQualities
        if (useReverseComplement)
            bq = bq.reversedArray()
        return bq
    }

    val baseQualityString: String get()
    {
        var bq = read.baseQualityString
        return if (useReverseComplement)
            bq.reversed()
        else bq
    }

    // this is the sequence that is the potential CDR3 sequence
    val anchorSequence: String get()
    {
        return readSequence.substring(Math.max(anchorOffsetStart, 0), Math.min(anchorOffsetEnd, readLength))
    }

    val anchorAA: String get()
    {
        return Codons.aminoAcidFromBases(anchorSequence)
    }
}

// this stores the combined V and J match
data class Cdr3ReadVJMatch(
    val vMatch: VJReadCandidate?,
    val jMatch: VJReadCandidate?,
    val sequence: String,
    val vSequenceOffset: Int = -1,
    val jSequenceOffset: Int = -1
)
{
    // return the non null match, if both are not null then return the
    // non supplementary one
    val mainMatch: VJReadCandidate get()
    {
        if (vMatch != null && jMatch != null)
            return if (vMatch.read.isSecondaryOrSupplementary) jMatch else vMatch
        return vMatch ?: jMatch!!
    }

    val vAnchorOffsetStart: Int get()
    {
        return if (vMatch != null) vMatch.anchorOffsetStart + vSequenceOffset else -1
    }

    val vAnchorOffsetEnd: Int get()
    {
        return if (vMatch != null) vMatch.anchorOffsetEnd + vSequenceOffset else -1
    }

    val jAnchorOffsetStart: Int get()
    {
        return if (jMatch != null) jMatch.anchorOffsetStart + jSequenceOffset else -1
    }

    val jAnchorOffsetEnd: Int get()
    {
        return if (jMatch != null) jMatch.anchorOffsetEnd + jSequenceOffset else -1
    }
}