package com.hartwig.hmftools.cider

import com.hartwig.hmftools.cider.CiderConstants.MATCHES_REF_KNOWN_CDR3_AA
import com.hartwig.hmftools.cider.annotation.AlignmentAnnotation
import com.hartwig.hmftools.cider.annotation.AlignmentStatus
import com.hartwig.hmftools.cider.primer.VdjPrimerMatch
import com.hartwig.hmftools.common.codon.Codons
import com.hartwig.hmftools.common.utils.IntPair
import htsjdk.samtools.SAMRecord
import org.apache.logging.log4j.LogManager
import java.util.*
import kotlin.collections.ArrayList
import kotlin.collections.HashMap

// some annotation that we want to print out
data class VdjAnnotation(val vdj: VDJSequence,
                         val filters: List<Filter>,
                         val vAlignedReads: Int,
                         val jAlignedReads: Int,
                         val vNonSplitReads: Int,
                         val jNonSplitReads: Int,
                         val cdr3SupportMin: Int,
                         val vSimilarityScore: Int?,
                         val jSimilarityScore: Int?,
                         val vPrimerMatchCount: Int,
                         val jPrimerMatchCount: Int,
                         var alignmentAnnotation: AlignmentAnnotation? = null)
{
    enum class Filter
    {
        PASS,
        MATCHES_REF,
        NO_V_ANCHOR,
        NO_J_ANCHOR,
        POOR_V_ANCHOR,
        POOR_J_ANCHOR,
        PARTIAL,
        NO_VDJ_ALIGNMENT,
        DUPLICATE,
        MIN_LENGTH,
        MAX_LENGTH,
        CDR3_DELETED,
        NO_HIGH_QUAL_SUPPORT
    }
    
    val passesFilter : Boolean get()
    {
        return filters.contains(Filter.PASS)
    }

    val locus : IgTcrLocus get()
    {
        return if (vdj.vAnchor != null)
            vdj.vAnchor.geneType.locus
        else
            vdj.jAnchor!!.geneType.locus
    }
}

// helper object to annotate the VDJ sequences that we found
class VdjAnnotator(private val adaptor: IVJReadLayoutAdaptor,
                   private val blosumSearcher: IAnchorBlosumSearcher)
{
    fun sortAndAnnotateVdjs(vdjSequences: List<VDJSequence>,
                            alignmentAnnotations: Collection<AlignmentAnnotation>,
                            primerMatches: List<VdjPrimerMatch>) : List<VdjAnnotation>
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
            val alignmentAnnotation: AlignmentAnnotation? = alignmentAnnotations.find { o -> o.vdjSequence === vdj }
            vdjAnnotations.add(annotateVdj(vdj, isDuplicate, alignmentAnnotation, primerMatches))
        }

        return vdjAnnotations
    }

    fun annotateVdj(vdj: VDJSequence, isDuplicate: Boolean, alignmentAnnotation: AlignmentAnnotation?,
                    primerMatches: List<VdjPrimerMatch>) : VdjAnnotation
    {
        val vAnchorByReadMatch: VJAnchorByReadMatch? = vdj.vAnchor as? VJAnchorByReadMatch
        val jAnchorByReadMatch: VJAnchorByReadMatch? = vdj.jAnchor as? VJAnchorByReadMatch
        val vAlignedReads: Int = vAnchorByReadMatch?.numReads ?: 0
        val jAlignedReads: Int = jAnchorByReadMatch?.numReads ?: 0
        val vNonSplitReads: Int = countNonSplitReads(vdj, VJ.V)
        val jNonSplitReads: Int = countNonSplitReads(vdj, VJ.J)
        val cdr3SupportMin: Int = calcCdr3SupportMin(vdj)

        val matchesRef: Boolean = vdjMatchesRef(vdj, vAlignedReads = vAlignedReads, jAlignedReads = jAlignedReads,
                                                vNonSplitReads = vNonSplitReads, jNonSplitReads = jNonSplitReads)

        val filterReasons = annotateFilterReasons(vdj, isDuplicate, matchesRef, cdr3SupportMin, alignmentAnnotation)

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
            vPrimerMatchCount = vPrimerMatchCount, jPrimerMatchCount = jPrimerMatchCount,
            alignmentAnnotation = alignmentAnnotation)
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
                vdj.sequence, vdj.jAnchor.geneType.pairedVjGeneTypes(),
                IAnchorBlosumSearcher.Mode.ALLOW_NEG_SIMILARITY,
                0,
                vdj.jAnchor.anchorBoundary)
        }
        else
        {
            if (vdj.vAnchor == null)
                return null

            anchorBlosumMatch = blosumSearcher.searchForAnchor(
                vdj.sequence, vdj.vAnchor.geneType.pairedVjGeneTypes(),
                IAnchorBlosumSearcher.Mode.ALLOW_NEG_SIMILARITY,
                vdj.vAnchor.anchorBoundary,
                vdj.length)
        }

        return anchorBlosumMatch
    }

    // TODO: probably simpler to just directly count the number of split reads instead of
    // trying to see if any split location straddle a boundary
    //
    // =====V anchor========----------------------------------=====J anchor========
    // |---------------------------------|        |-------------------------------|
    //        V anchor + 30 bases                      30 bases + J anchor
    // We want to count reads that has an aligned block that span the region V anchor + 30 bases after
    // or J anchor + 30 bases before
    //
    // This would tell us if there is an rearrangement that occurred.
    fun countNonSplitReads(vdj: VDJSequence, vj: VJ) : Int
    {
        val alignedPos = vdj.layout.alignedPosition - vdj.layoutSliceStart

        if (vj == VJ.V && vdj.vAnchor == null)
            return 0

        if (vj == VJ.J && vdj.jAnchor == null)
            return 0

        // we want to work out if any read block spans the boundary region
        // which is the boundary +- 30 bases
        var boundaryRegionStart: Int
        var boundaryRegionEnd: Int

        if (vj == VJ.V)
        {
            boundaryRegionStart = vdj.vAnchorBoundary!!
            boundaryRegionEnd = vdj.vAnchorBoundary!!
        }
        else
        {
            boundaryRegionStart = vdj.jAnchorBoundary!!
            boundaryRegionEnd = vdj.jAnchorBoundary!!
        }

        boundaryRegionStart -= CiderConstants.MIN_NON_SPLIT_READ_STRADDLE_LENGTH
        boundaryRegionEnd += CiderConstants.MIN_NON_SPLIT_READ_STRADDLE_LENGTH

        // if the VDJ sequence only get a partial anchor sequence, any align block that span the whole way is counted
        // as non split read
        // i.e. if the V anchor has only 10 bases, then any align block that spans that 10 bases + 30 bases after V boundary
        // is counted as non split read
        boundaryRegionStart = Math.max(boundaryRegionStart, 0)
        boundaryRegionEnd = Math.min(boundaryRegionEnd, vdj.length)

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

            // we are trying to find if there is an aligned block that span the following region
            val readSliceBoundaryRegionStart = boundaryRegionStart - readPosWithinVdj
            val readSliceBoundaryRegionEnd = boundaryRegionEnd - readPosWithinVdj

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
                    if (readSliceBoundaryRegionStart >= alignRangeInReadSlice.left &&
                        readSliceBoundaryRegionEnd <= alignRangeInReadSlice.right
                    )
                    {
                        ++nonSplitReadCount

                        // this read straddles an anchor boundary
                        sLogger.trace("read({}) cigar({}) revcomp({}), straddles {} boundary, align offset({}:{}), boundary region({}:{})",
                            samRecord, samRecord.cigarString, layoutReadSlice.reverseComplement, vj,
                            alignRangeInReadSlice.left, alignRangeInReadSlice.right,
                            readSliceBoundaryRegionStart, readSliceBoundaryRegionEnd)
                    }
                }
            }
        }
        return nonSplitReadCount
    }

    fun vdjMatchesRef(vdj: VDJSequence, vAlignedReads: Int, jAlignedReads: Int,
                      vNonSplitReads: Int, jNonSplitReads: Int) : Boolean
    {
        val maxNonSplitReads = Math.min(MAX_NONSPLIT_READS, Math.max(vdj.numReads / 2, 1))
        return ((vAlignedReads == 0 || jAlignedReads == 0) &&
            (jNonSplitReads + vNonSplitReads) >= maxNonSplitReads)
                || vdjInKnownMatchesRefList(vdj)
    }

    fun vdjMatchesRef(vdj: VDJSequence) : Boolean
    {
        val vAnchorByReadMatch: VJAnchorByReadMatch? = vdj.vAnchor as? VJAnchorByReadMatch
        val jAnchorByReadMatch: VJAnchorByReadMatch? = vdj.jAnchor as? VJAnchorByReadMatch
        val vAlignedReads: Int = vAnchorByReadMatch?.numReads ?: 0
        val jAlignedReads: Int = jAnchorByReadMatch?.numReads ?: 0
        val vNonSplitReads: Int = countNonSplitReads(vdj, VJ.V)
        val jNonSplitReads: Int = countNonSplitReads(vdj, VJ.J)
        return vdjMatchesRef(vdj, vAlignedReads = vAlignedReads, jAlignedReads = jAlignedReads,
            vNonSplitReads = vNonSplitReads, jNonSplitReads = jNonSplitReads)
    }

    fun vdjInKnownMatchesRefList(vdj: VDJSequence): Boolean
    {
        val cdr3AA = Codons.aminoAcidFromBases(vdj.cdr3Sequence)
        return MATCHES_REF_KNOWN_CDR3_AA.any { seq -> cdr3AA.matches(seq) }
    }

    companion object
    {
        const val CDR3_FILTER_AA_MIN_LENGTH = 5
        const val CDR3_FILTER_AA_MAX_LENGTH = 40
        const val MAX_NONSPLIT_READS = 2

        private val sLogger = LogManager.getLogger(VdjAnnotator::class.java)

        fun calcCdr3SupportMin(vdj: VDJSequence) : Int
        {
            return vdj.supportCounts.minOrNull() ?: 0
        }

        fun calcAnchorSimilarity(vdj: VDJSequence, anchor: VJAnchor) : Int
        {
            val seq: String
            val templateAnchorSeq: String = anchor.templateAnchorSeq

            if (anchor.vj == VJ.V)
            {
                seq = vdj.vAnchorSequence!!
            }
            else
            {
                seq = vdj.jAnchorSequence!!
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

        fun annotateFilterReasons(vdj: VDJSequence, isDuplicate: Boolean, matchesRef: Boolean, cdr3SupportMin: Int,
                                  alignmentAnnotation: AlignmentAnnotation?): List<VdjAnnotation.Filter>
        {
            val filters = ArrayList<VdjAnnotation.Filter>()

            if (matchesRef || alignmentAnnotation?.alignmentStatus == AlignmentStatus.NO_REARRANGEMENT)
            {
                filters.add(VdjAnnotation.Filter.MATCHES_REF)
            }
            if (alignmentAnnotation == null)
            {
                // following filters are only applied if we do not use alignment
                if (vdj.vAnchor == null)
                {
                    filters.add(VdjAnnotation.Filter.NO_V_ANCHOR)
                }
                if (vdj.jAnchor == null)
                {
                    filters.add(VdjAnnotation.Filter.NO_J_ANCHOR)
                }
                if ((vdj.vAnchor is VJAnchorByBlosum) && (vdj.vAnchor.similarityScore < 0))
                {
                    filters.add(VdjAnnotation.Filter.POOR_V_ANCHOR)
                }
                if ((vdj.jAnchor is VJAnchorByBlosum) && (vdj.jAnchor.similarityScore < 0))
                {
                    filters.add(VdjAnnotation.Filter.POOR_J_ANCHOR)
                }
            }
            else
            {
                if (alignmentAnnotation.alignmentStatus == AlignmentStatus.V_D ||
                    alignmentAnnotation.alignmentStatus == AlignmentStatus.D_J ||
                    alignmentAnnotation.alignmentStatus == AlignmentStatus.V_ONLY ||
                    alignmentAnnotation.alignmentStatus == AlignmentStatus.J_ONLY)
                {
                    filters.add(VdjAnnotation.Filter.PARTIAL)
                }
                else if (alignmentAnnotation.alignmentStatus == AlignmentStatus.NO_VDJ_ALIGNMENT)
                {
                    filters.add(VdjAnnotation.Filter.NO_VDJ_ALIGNMENT)
                }
            }
            if (isDuplicate)
            {
                filters.add(VdjAnnotation.Filter.DUPLICATE)
            }
            if (vdj.cdr3Sequence.length < CDR3_FILTER_AA_MIN_LENGTH * 3)
            {
                filters.add(VdjAnnotation.Filter.MIN_LENGTH)
            }
            if (vdj.cdr3Sequence.length > CDR3_FILTER_AA_MAX_LENGTH * 3)
            {
                filters.add(VdjAnnotation.Filter.MAX_LENGTH)
            }
            if (vdj.vAnchor != null && vdj.jAnchor != null &&
                vdj.vAnchorBoundary!! > vdj.jAnchorBoundary!!)
            {
                filters.add(VdjAnnotation.Filter.CDR3_DELETED)
            }
            if (cdr3SupportMin <= 0)
            {
                filters.add(VdjAnnotation.Filter.NO_HIGH_QUAL_SUPPORT)
            }
            if (filters.isEmpty())
            {
                filters.add(VdjAnnotation.Filter.PASS)
            }
            return filters
        }
    }
}
