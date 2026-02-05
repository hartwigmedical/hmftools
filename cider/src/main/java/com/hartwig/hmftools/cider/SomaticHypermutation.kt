package com.hartwig.hmftools.cider

import com.hartwig.hmftools.cider.annotation.ImgtSequenceFile
import com.hartwig.hmftools.cider.genes.VJ
import com.hartwig.hmftools.common.genome.region.Strand
import htsjdk.samtools.util.SequenceUtil.reverseComplement

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
    val indelBases: Int
)


// TODO: unit test
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
    val layoutSeqStranded = if (isForward) layoutSeq else reverseComplement(layoutSeq)
    val layoutStart = if (isForward) queryRange.start else layoutSeq.length - 1 - queryRange.endInclusive
    val layoutSeqBounds = if (isV)
        (if (isForward) 0 until layoutAnchorBoundary else (layoutSeq.length - layoutAnchorBoundary) until layoutSeq.length)
    else (if (isForward) layoutAnchorBoundary until layoutSeq.length else 0 until layoutSeq.length - layoutAnchorBoundary)

    val imgtSeq = imgtSequence.sequenceWithRef
    val imgtAlignStart = alignment.refRange.start
    val imgtSeqBounds = imgtSequence.imgtRange

    var comparedBases = 0
    var matches = 0
    var indelBases = 0
    var layoutIndex = layoutStart
    var imgtIndex = imgtAlignStart
    var firstComparedLayoutIndex: Int? = null
    var lastComparedLayoutIndex = 0
    var firstComparedImgtIndex: Int? = null
    var lastComparedImgtIndex = 0

    for (element in alignment.cigar)
    {
        val op = element.operator
        if (op.isAlignment)
        {
            for (i in 0 until element.length)
            {
                if (layoutSeqBounds.contains(layoutIndex) && imgtSeqBounds.contains(imgtIndex))
                {
                    comparedBases++
                    if (firstComparedLayoutIndex == null)
                    {
                        firstComparedLayoutIndex = layoutIndex
                    }
                    lastComparedLayoutIndex = layoutIndex
                    if (firstComparedImgtIndex == null)
                    {
                        firstComparedImgtIndex = imgtIndex
                    }
                    lastComparedImgtIndex = imgtIndex

                    if (layoutSeqStranded[layoutIndex] == imgtSeq[imgtIndex])
                    {
                        matches++
                    }
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

    return if (comparedBases == 0) null
    else ShmGeneComparison(
        seqLength = lastComparedLayoutIndex - firstComparedLayoutIndex!! + 1,
        imgtLength = lastComparedImgtIndex - firstComparedImgtIndex!! + 1,
        // For now, we are excluding indels in the %identity calculation because we're not sure how they should be counted, and
        // there are usually few of them
        // TODO? handle indels better
        pctIdentity = (100.0 * matches) / comparedBases,
        indelBases = indelBases)
}
