package com.hartwig.hmftools.cider

import com.hartwig.hmftools.cider.primer.VdjPrimerMatch
import com.hartwig.hmftools.common.utils.IntPair
import htsjdk.samtools.SAMRecord
import org.apache.logging.log4j.LogManager
import java.util.*
import kotlin.collections.ArrayList
import kotlin.collections.HashMap

// some annotation that we want to print out
data class VdjAnnotation(val vdj: VDJSequence,
                         val filters: List<String>,
                         val vAlignedReads: Int,
                         val jAlignedReads: Int,
                         val vNonSplitReads: Int,
                         val jNonSplitReads: Int,
                         val cdr3SupportMin: Int,
                         val vSimilarityScore: Int?,
                         val jSimilarityScore: Int?,
                         val vPrimerMatchCount: Int,
                         val jPrimerMatchCount: Int)

// helper object to annotate the VDJ sequences that we found
class VdjAnnotator(private val adaptor: VJReadLayoutBuilder,
                   private val blosumSearcher: AnchorBlosumSearcher)
{
    fun sortAndAnnotateVdjs(vdjSequences: List<VDJSequence>, primerMatches: List<VdjPrimerMatch>) : List<VdjAnnotation>
    {
        val sortedVdj = vdjSequences.sortedWith(
            Collections.reverseOrder(
                Comparator.comparingInt({ vdj: VDJSequence -> calcCdr3SupportMin(vdj) })
                    .thenComparingInt({ vdj: VDJSequence -> vdj.numReads }) // handle the highest quality ones first
                    .thenComparingInt({ vdj: VDJSequence -> vdj.length })
                    .thenComparing({ vdj: VDJSequence -> vdj.cdr3Sequence })
            ))

        // we create a set of all VDJs that are NOT duplicates
        // keyed by the cdr3Sequence
        val notDuplicateVdjs = HashMap<String, VDJSequence>()

        for (vdjSeq in sortedVdj)
        {
            val cdr3 = vdjSeq.cdr3Sequence

            val vdjWithMoreSupport : VDJSequence? = notDuplicateVdjs[cdr3]

            if (vdjWithMoreSupport == null)
            {
                notDuplicateVdjs[cdr3] = vdjSeq
            }
            else
            {
                // just make sure it has more support, this should be guaranteed by the
                // previous sort
                require(calcCdr3SupportMin(vdjWithMoreSupport) >= calcCdr3SupportMin(vdjSeq))
            }
        }

        val vdjAnnotations = ArrayList<VdjAnnotation>()

        for (vdj in sortedVdj)
        {
            val vdjWithMoreSupport : VDJSequence? = notDuplicateVdjs[vdj.cdr3Sequence]
            requireNotNull(vdjWithMoreSupport)
            val isDuplicate: Boolean = vdjWithMoreSupport !== vdj
            vdjAnnotations.add(annotateVdj(vdj, isDuplicate, primerMatches))
        }

        return vdjAnnotations
    }

    fun annotateVdj(vdj: VDJSequence, isDuplicate: Boolean, primerMatches: List<VdjPrimerMatch>) : VdjAnnotation
    {
        val vAnchorByReadMatch: VJAnchorByReadMatch? = vdj.vAnchor as? VJAnchorByReadMatch
        val jAnchorByReadMatch: VJAnchorByReadMatch? = vdj.jAnchor as? VJAnchorByReadMatch
        val vAlignedReads: Int = vAnchorByReadMatch?.numReads ?: 0
        val jAlignedReads: Int = jAnchorByReadMatch?.numReads ?: 0
        val vNonSplitReads: Int = countNonSplitReads(vdj, VJ.V)
        val jNonSplitReads: Int = countNonSplitReads(vdj, VJ.J)
        val cdr3SupportMin: Int = calcCdr3SupportMin(vdj)
        val filterReasons = annotateFilterReasons(vdj, isDuplicate, vAlignedReads, jAlignedReads, vNonSplitReads, jNonSplitReads)

        val vSimilarityScore: Int? = if (vdj.vAnchor != null) calcAnchorSimilarity(vdj, vdj.vAnchor) else null
        val jSimilarityScore: Int? = if (vdj.jAnchor != null) calcAnchorSimilarity(vdj, vdj.jAnchor) else null

        // also count how many primer matches
        // look through all the primer matches and see which ones has this vdj in it
        val vPrimerMatchCount = primerMatches.filter({ o -> o.vdj === vdj && o.primer.vj == VJ.V}).size
        val jPrimerMatchCount = primerMatches.filter({ o -> o.vdj === vdj && o.primer.vj == VJ.J}).size

        return VdjAnnotation(vdj, filterReasons,
            vAlignedReads = vAlignedReads, jAlignedReads = jAlignedReads,
            vNonSplitReads = vNonSplitReads, jNonSplitReads = jNonSplitReads,
            cdr3SupportMin = cdr3SupportMin,
            vSimilarityScore = vSimilarityScore, jSimilarityScore = jSimilarityScore,
            vPrimerMatchCount = vPrimerMatchCount, jPrimerMatchCount = jPrimerMatchCount)
    }

    fun findAnchorByBlosum(vdj: VDJSequence, vj: VJ) : AnchorBlosumMatch?
    {
        // see if we can find an anchor in this sequence
        val anchorBlosumMatch: AnchorBlosumMatch?
        if (vj == VJ.V)
        {
            if (vdj.jAnchor == null)
                return null

            anchorBlosumMatch = blosumSearcher.searchForAnchor(
                vdj.sequence, vdj.jAnchor.geneType.pairedVjGeneType(),
                IAnchorBlosumSearcher.Mode.ALLOW_NEG_SIMILARITY,
                0,
                vdj.jAnchor.anchorBoundary)
        }
        else
        {
            if (vdj.vAnchor == null)
                return null

            anchorBlosumMatch = blosumSearcher.searchForAnchor(
                vdj.sequence, vdj.vAnchor.geneType.pairedVjGeneType(),
                IAnchorBlosumSearcher.Mode.ALLOW_NEG_SIMILARITY,
                vdj.vAnchor.anchorBoundary,
                vdj.length)
        }

        return anchorBlosumMatch
    }

    // for now we want to just log where the mappings are
    fun countNonSplitReads(vdj: VDJSequence, vj: VJ) : Int
    {
        val alignedPos = vdj.layout.alignedPosition - vdj.layoutSliceStart

        if (vj == VJ.V && vdj.vAnchor == null)
            return 0

        if (vj == VJ.J && vdj.jAnchor == null)
            return 0

        val boundaryPos: Int = if (vj == VJ.V)
        {
            vdj.vAnchorBoundary!!
        }
        else
        {
            vdj.jAnchorBoundary!!
        }
        var nonSplitReadCount = 0

        for (read in vdj.layout.reads)
        {
            val samRecord: SAMRecord = adaptor.toReadCandidate(read).read
            val layoutReadSlice: ReadSlice = adaptor.toLayoutReadSlice(read)

            // AGATCTGAG-GACACGGCCGTGTATTACTGT-GCGAGAGACACAGTGTGAAAACCCACATCCTGAGAGTGTCAGAAACCCTGAGGGA
            //           |___________________|
            //                 V anchor
            //                  =========================> read slice
            //                                        |  <-- aligned position
            //                  |------------|           <-- this is the value we want
            val readPosWithinVdj = alignedPos - read.alignedPosition
            val readSliceAnchorBoundary = boundaryPos - readPosWithinVdj

            // work out where in the VDJ sequence is this read mapped
            if (!samRecord.readUnmappedFlag)
            {
                for (alignBlock in CiderUtils.getAdjustedAlignmentBlocks(samRecord.cigar))
                {
                    // now get those positions in terms of read slice
                    val alignRangeInReadSlice: IntPair = layoutReadSlice.readRangeToSliceRange(
                        alignBlock.readStart - 1,
                        alignBlock.readStart - 1 + alignBlock.length)

                    require(alignRangeInReadSlice.left < alignRangeInReadSlice.right)

                    // to count as a non split read, it needs to be away from the boundary
                    if (readSliceAnchorBoundary >= alignRangeInReadSlice.left + CiderConstants.MIN_NON_SPLIT_READ_STRADDLE_LENGTH &&
                        readSliceAnchorBoundary <= alignRangeInReadSlice.right - CiderConstants.MIN_NON_SPLIT_READ_STRADDLE_LENGTH
                    )
                    {
                        ++nonSplitReadCount

                        // this read straddles a v anchor boundary
                        sLogger.trace("read({}) cigar({}) revcomp({}), straddles {} boundary, align offset({}:{}), boundary offset({})",
                            samRecord, samRecord.cigarString, layoutReadSlice.reverseComplement, vj,
                            alignRangeInReadSlice.left, alignRangeInReadSlice.right, readSliceAnchorBoundary)
                    }
                }
            }
        }
        return nonSplitReadCount
    }

    companion object
    {
        const val CDR3_FILTER_AA_MIN_LENGTH = 5
        const val CDR3_FILTER_AA_MAX_LENGTH = 40
        const val MAX_NONSPLIT_READS = 2

        private val sLogger = LogManager.getLogger(VdjAnnotator::class.java)

        fun calcCdr3SupportMin(vdj: VDJSequence) : Int
        {
            val supportCounts = vdj.supportCounts
            var supportSliceStart: Int = vdj.cdr3Start
            var supportSliceEnd: Int = vdj.cdr3End

            // if V or J anchor is missing we want to limit them to just 60 bases
            if (vdj.vAnchor == null)
            {
                supportSliceStart = Math.max(supportSliceEnd - CiderConstants.PARTIAL_VDJ_UNANCHORED_LENGTH_BASES, 0)
            }
            else if (vdj.jAnchor == null)
            {
                supportSliceEnd = Math.min(supportSliceStart + CiderConstants.PARTIAL_VDJ_UNANCHORED_LENGTH_BASES, supportCounts.size)
            }
            return supportCounts.slice(supportSliceStart until supportSliceEnd).minOrNull() ?: 0
        }

        fun calcAnchorSimilarity(vdj: VDJSequence, anchor: VJAnchor) : Int
        {
            val seq: String
            val templateAnchorSeq: String = anchor.templateAnchorSeq

            if (anchor.vj == VJ.V)
            {
                seq = vdj.vAnchorSequence
            }
            else
            {
                seq = vdj.jAnchorSequence
            }

            try
            {
                return BlosumSimilarityCalc.calcSimilarityScore(anchor.vj, templateAnchorSeq, seq)
            }
            catch (e: IllegalArgumentException)
            {
                throw IllegalArgumentException("cannot calc similarity score: ${templateAnchorSeq} and ${seq}")
            }
        }

        fun annotateFilterReasons(vdj: VDJSequence, isDuplicate: Boolean,
                         vAlignedReads: Int, jAlignedReads: Int,
                        vNonSplitReads: Int, jNonSplitReads: Int): List<String>
        {
            val filters = ArrayList<String>()

            val maxNonSplitReads = Math.min(MAX_NONSPLIT_READS, Math.max(vdj.numReads / 2, 1))
            if ((vAlignedReads == 0 || jAlignedReads == 0) &&
                (jNonSplitReads + vNonSplitReads) >= maxNonSplitReads)
            {
                filters.add("MATCHES_REF")
            }
            if (vdj.vAnchor == null)
            {
                filters.add("NO_V_ANCHOR")
            }
            if (vdj.jAnchor == null)
            {
                filters.add("NO_J_ANCHOR")
            }
            if ((vdj.vAnchor is VJAnchorByBlosum) && (vdj.vAnchor.similarityScore < 0))
            {
                filters.add("POOR_V_ANCHOR")
            }
            if ((vdj.jAnchor is VJAnchorByBlosum) && (vdj.jAnchor.similarityScore < 0))
            {
                filters.add("POOR_J_ANCHOR")
            }
            if (isDuplicate)
            {
                filters.add("DUPLICATE")
            }
            if (vdj.cdr3Sequence.length < CDR3_FILTER_AA_MIN_LENGTH * 3)
            {
                filters.add("MIN_LENGTH")
            }
            if (vdj.cdr3Sequence.length > CDR3_FILTER_AA_MAX_LENGTH * 3)
            {
                filters.add("MAX_LENGTH")
            }
            if (filters.isEmpty())
            {
                filters.add("PASS")
            }
            return filters
        }
    }
}
