package com.hartwig.hmftools.cider

import com.hartwig.hmftools.cider.annotation.ImgtSequenceFile
import com.hartwig.hmftools.cider.genes.VJ
import com.hartwig.hmftools.common.genome.region.Strand
import htsjdk.samtools.util.SequenceUtil.reverseComplement
import kotlin.math.max
import kotlin.math.min

// Determination of "somatic hypermutation status" from V-region sequence, which is a prognostic indicator for chronic lymphocytic leukemia.
// https://pmc.ncbi.nlm.nih.gov/articles/PMC7248390/


enum class SomaticHypermutationStatus {
    UNMUTATED,
    MUTATED_BORDERLINE,
    MUTATED;

    companion object
    {
        // Determine the status from the % identity of the V-region to the IMGT reference sequence.
        fun determineFromVGeneIdentity(pctIdentity: Double): SomaticHypermutationStatus
        {
            return if (pctIdentity >= 98) UNMUTATED
            else if (pctIdentity >= 97) MUTATED_BORDERLINE
            else MUTATED
        }
    }
}


data class ShmGeneComparison(
    val seqLength: Int,         // Length of the VDJ subsequence used in comparison.
    val imgtLength: Int,        // Length of the IMGT sequence used in comparison.
    val pctIdentity: Double,    // Range [0, 100]. Excludes indel bases (for now, at least).
    val indelBases: Int,        // Number of bases in indels in the alignment.
    val clipBases: Int          // Number of bases in the IMGT sequence which could've been aligned but were clipped.
)
{
    init
    {
        require(seqLength >= 0)
        require(imgtLength >= 0)
        require(pctIdentity >= 0.0 && pctIdentity <= 100.0)
        require(indelBases >= 0)
        require(clipBases >= 0)
    }
}


fun compareVJRegionToImgt(
    layoutSeq: String, type: VJ, layoutAnchorBoundary: Int, imgtSequence: ImgtSequenceFile.Sequence, queryRange: IntRange,
    alignment: Alignment
): ShmGeneComparison?
{
    // Calculate percentage identity of the V region between the sample and IMGT sequence, which is a heuristic for the degree of somatic
    // hypermutation.
    // The V region is from Cys104 (last amino acid of anchor) upstream to the start of the second exon.
    // It seems like the start of the sequence in the IMGT resource is the start of the second exon.

    //                           V side      || anc ||    CDR3    || anc ||  J side
    // layout:          LLLLLLLLLLLLLLLLLLLLLLAAAAAAACCCCCCCCCCCCCCAAAAAAALLLLLLLLLLLLLLLLLL
    // alignment query:       |--------------------------------------------------------|
    // IMGT sequence:     RRRRRRRRRRIIIIIIIIIIIIIIIIIIIIIIIII
    //                    |  ref   ||         IMGT          |
    // alignment:                |-----------------------|
    // compare:                     |---------------|

    // We do a similar logic for the J side for completeness, even though that is not currently associated with a clinical marker.

    val isV = type == VJ.V
    val isForward = alignment.strand == Strand.FORWARD
    val layoutStart = queryRange.start
    val layoutEnd = queryRange.endInclusive + 1
    val layoutSeqBounds = if (isV) 0 until layoutAnchorBoundary else (layoutAnchorBoundary until layoutSeq.length)

    // It's easiest to conceptualise if we apply the strand to the ref (IMGT) rather than query (layout), even though usually alignment is
    // relative to the ref forward strand.
    // That way we can keep the layout in transcription order.
    val imgtSeq = if (isForward) imgtSequence.sequenceWithRef else reverseComplement(imgtSequence.sequenceWithRef)
    val imgtAlignStart = if (isForward) alignment.refRange.start else imgtSeq.length - (alignment.refRange.endInclusive + 1)
    val imgtAlignEnd = if (isForward) alignment.refRange.endInclusive + 1 else imgtSeq.length - alignment.refRange.start
    val imgtSeqBounds = if (isForward) imgtSequence.imgtRange
        else imgtSeq.length - (imgtSequence.imgtRange.endInclusive + 1) until imgtSeq.length - imgtSequence.imgtRange.start

    // CIGAR is always relative to the ref forward strand, so need to reverse it if we're using the ref reverse strand.
    val cigar = if (isForward) alignment.cigar else alignment.cigar.reversed()

    var comparedBases = 0
    var matches = 0
    var indelBases = 0
    var layoutIndex = layoutStart
    var imgtIndex = imgtAlignStart
    var firstAlignedLayoutIndex: Int? = null
    var firstComparedLayoutIndex: Int? = null
    var firstComparedImgtIndex: Int? = null
    var lastAlignedLayoutIndex = 0
    var lastComparedLayoutIndex = 0
    var lastComparedImgtIndex = 0

    for (element in cigar)
    {
        val op = element.operator
        if (op.isAlignment)
        {
            firstAlignedLayoutIndex = firstAlignedLayoutIndex ?: layoutIndex
            lastAlignedLayoutIndex = layoutIndex + element.length - 1
            for (i in 0 until element.length)
            {
                if (layoutSeqBounds.contains(layoutIndex) && imgtSeqBounds.contains(imgtIndex))
                {
                    comparedBases++
                    matches += if (layoutSeq[layoutIndex] == imgtSeq[imgtIndex]) 1 else 0
                    firstComparedLayoutIndex = firstComparedLayoutIndex ?: layoutIndex
                    firstComparedImgtIndex = firstComparedImgtIndex ?: imgtIndex
                    lastComparedLayoutIndex = layoutIndex
                    lastComparedImgtIndex = imgtIndex
                }
                layoutIndex += 1
                imgtIndex += 1
            }
        }
        else if (op.consumesReadBases() || op.isClipping)
        {
            layoutIndex += element.length
        }
        else if (op.consumesReferenceBases())
        {
            imgtIndex += element.length
        }
        if (op.isIndel)
        {
            indelBases += element.length
        }
    }

    assert(layoutIndex == layoutEnd)
    assert(imgtIndex == imgtAlignEnd)

    return if (comparedBases == 0) null
    else ShmGeneComparison(
        seqLength = lastComparedLayoutIndex - firstComparedLayoutIndex!! + 1,
        imgtLength = lastComparedImgtIndex - firstComparedImgtIndex!! + 1,
        // For now, we are excluding indels in the %identity calculation because we're not sure how they should be counted, and
        // there are usually few of them.
        pctIdentity = 100.0 * matches / comparedBases,
        indelBases = indelBases,
        clipBases = max(0,
            if (isV) min(firstAlignedLayoutIndex!! - layoutStart, imgtAlignStart - imgtSeqBounds.start)
                    else min(layoutEnd - lastAlignedLayoutIndex, imgtSeqBounds.endInclusive + 1 - imgtAlignEnd)))
}
