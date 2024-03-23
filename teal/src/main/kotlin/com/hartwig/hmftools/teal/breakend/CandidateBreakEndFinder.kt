package com.hartwig.hmftools.teal.breakend

import com.hartwig.hmftools.common.genome.chromosome.ContigComparator
import com.hartwig.hmftools.common.bam.CigarUtils
import htsjdk.samtools.SAMRecord
import com.hartwig.hmftools.common.region.ExcludedRegions
import com.hartwig.hmftools.common.region.ChrBaseRegion
import com.hartwig.hmftools.teal.ReadGroup
import com.hartwig.hmftools.teal.TealUtils
import com.hartwig.hmftools.teal.util.TelomereMatcher
import org.apache.logging.log4j.LogManager
import java.util.Comparator
import kotlin.math.abs

const val BREAKEND_MERGE_DISTANCE = 2

// go through the bam file and find any potential
// split read sites
class CandidateBreakEndFinder(
    config: BreakEndParams,
    extraCandidateBreakEnds: Collection<TelomericBreakEnd> = ArrayList())
{
    private val logger = LogManager.getLogger(javaClass)
    
    data class CandidateBreakEnd(
        val key: TelomericBreakEnd,
        val splitReads: MutableList<SAMRecord> = ArrayList())
    {
        constructor(type: TelomericBreakEndType,
                    chromosome: String,
                    position: Int
        ) : this(TelomericBreakEnd(type, chromosome, position))

        val type: TelomericBreakEndType get() = key.type
        val chromosome: String get() = key.chromosome
        val position: Int get() = key.position
    }

    private val mConfig = config

    /*
    // we let some bases at the edge to be non telomeric
    private val mSplitTelomereMatchThreshold = splitTelomereMatchThreshold

    // I find using 0.8 can still get all the break ends gripss found.
    // using 0.7 and we start excluding some
    private val mAlignedTelomereMatchThreshold = alignedTelomereMatchThreshold
    private val mFuzzyMatchDistance = fuzzyMatchDistance
    private val mMaxAlignedPolyBaseThreshold = maxAlignedPolyBaseThreshold
    private val mExcludedGenomeRegion = excludedGenomeRegion
     */

    private var mBoundaryZone = 1
    private val mMinTelomereMatchLength = 12

    private val mCandidateBreakEnds: MutableList<CandidateBreakEnd> = ArrayList()
    val potentialBreakEnds: List<TelomericBreakEnd> get() = mCandidateBreakEnds.map({ o -> o.key })

    init
    {
        // we need to add the candidate breakends
        extraCandidateBreakEnds.forEach({ k -> mCandidateBreakEnds.add(CandidateBreakEnd(k)) })
    }

    fun setBoundaryZone(bz: Int)
    {
        mBoundaryZone = bz
    }

    private fun matchesGTelomere(seq: String, telomereMatchThreshold: Double): Boolean
    {
        // to speed it up we do a quick count first
        return TelomereMatcher.matchesGTelomere(seq, telomereMatchThreshold, mMinTelomereMatchLength)
    }

    private fun matchesCTelomere(seq: String, telomereMatchThreshold: Double): Boolean
    {
        return TelomereMatcher.matchesCTelomere(seq, telomereMatchThreshold, mMinTelomereMatchLength)
    }

    fun processReadGroups(readGroups: Collection<ReadGroup>)
    {
        for (rg in readGroups)
        {
            for (r in rg.allReads)
            {
                // for now we just process the records without any relation
                // the the group. Soon we might want to change the way we
                // do it
                processRead(r)
            }
        }
        consolidatePotentialBreakEnds()
    }

    // we want to collect evidence to support position of a break end
    fun processRead(r: SAMRecord)
    {
        if (r.readUnmappedFlag)
            return

        val tbe = findCandidateBreakEnd(r)
        if (tbe != null)
        {
            mCandidateBreakEnds.add(tbe)
            //TE_LOGGER.trace("record: {} breakend: {} cigar: {} readBases: {}", r, tbe, r.getCigarString(), r.getReadString());
        }
    }

    // look at soft clip site, and one part is mapped to a genome, and the
    // other part is telomeric
    // A right soft clip looks like 73S78M
    // A left soft clip looks like 78M32S
    // we want to work out whether one side is telomeric and the other side is not
    fun findCandidateBreakEnd(r: SAMRecord): CandidateBreakEnd?
    {
        // TE_LOGGER.trace("record: {}", r);
        if (isInExcludedBaseRegion(mConfig, r.referenceName, r.alignmentStart, r.alignmentEnd))
        {
            return null
        }
        val readString = r.readString
        val leftClip = CigarUtils.leftSoftClipLength(r)
        val rightClip = CigarUtils.rightSoftClipLength(r)

        // now we work out if the clipped part is telomere
        val leftClipBases = readString.substring(0, leftClip)
        val rightClipBases = readString.substring(readString.length - rightClip)
        val alignedBases = readString.substring(leftClip, readString.length - rightClip)

        // check poly A / T threshold
        val polyBaseRatio = maxOf(alignedBases.count({c -> c == 'A'}),
                                  alignedBases.count({c -> c == 'C'}),
                                  alignedBases.count({c -> c == 'G'}),
                                  alignedBases.count({c -> c == 'T'})).toDouble() / alignedBases.length.toDouble()

        if (polyBaseRatio > mConfig.maxAlignedPolyBaseThreshold)
        {
            return null
        }

        assert(leftClipBases.length == leftClip)
        assert(rightClipBases.length == rightClip)
        assert(leftClipBases.length + rightClipBases.length + alignedBases.length == readString.length)

        if (alignedBases.length <= mMinTelomereMatchLength)
        {
            // if the aligned section is too short to tell whether it is
            // even telomeric or not
            return null
        }

        if (TealUtils.isPolyGC(alignedBases))
        {
            // aligned part must not be poly GC
            return null
        }

        logger.trace("checking record: {}", r)

        var tbe: CandidateBreakEnd? = null
        if (leftClipBases.length >= mMinTelomereMatchLength + mBoundaryZone)
        {
            // we remove boundary zone
            val modLeftClipBases = leftClipBases.substring(0, leftClipBases.length - mBoundaryZone)
            if (matchesGTelomere(modLeftClipBases, mConfig.telomereMatchThreshold))
            {
                tbe = CandidateBreakEnd(TelomericBreakEndType.LEFT_G_TELOMERIC, r.referenceName, r.alignmentStart)
            }
            else if (matchesCTelomere(modLeftClipBases, mConfig.telomereMatchThreshold))
            {
                tbe = CandidateBreakEnd(TelomericBreakEndType.LEFT_C_TELOMERIC, r.referenceName, r.alignmentStart)
            }
        }
        if (tbe == null && rightClipBases.length >= mMinTelomereMatchLength + mBoundaryZone)
        {
            // we remove boundary zone
            val modRightClipBases = rightClipBases.substring(mBoundaryZone)
            if (matchesGTelomere(modRightClipBases, mConfig.telomereMatchThreshold))
            {
                tbe = CandidateBreakEnd(TelomericBreakEndType.RIGHT_G_TELOMERIC, r.referenceName, r.alignmentEnd)
            } else if (matchesCTelomere(modRightClipBases, mConfig.telomereMatchThreshold))
            {
                tbe = CandidateBreakEnd(TelomericBreakEndType.RIGHT_C_TELOMERIC, r.referenceName, r.alignmentEnd)
            }
        }

        if (tbe == null)
        {
            return null
        }

        // this could be potential telomeric region
        // next thing we need to check that the aligned region is not the same type of telomere
        when (tbe.key.type)
        {
            TelomericBreakEndType.LEFT_C_TELOMERIC, TelomericBreakEndType.RIGHT_C_TELOMERIC ->
            {
                if (matchesCTelomere(alignedBases, mConfig.alignedSegmentTelomereMatchThreshold))
                {
                    // aligned bases is also  C telomeric, so this is no longer interesting
                    tbe = null
                }
            }
            TelomericBreakEndType.LEFT_G_TELOMERIC, TelomericBreakEndType.RIGHT_G_TELOMERIC ->
            {
                if (matchesGTelomere(alignedBases, mConfig.alignedSegmentTelomereMatchThreshold))
                {
                    // aligned bases is also  C telomeric, so this is no longer interesting
                    tbe = null
                }
            }
        }

        if (tbe != null)
        {
            // add this to list of split reads
            tbe.splitReads.add(r)
            logger.trace(
                "record({}) -v strand({}) breakend({}) cigar({}) leftClip({}) rightClip({}) aligned({})",
                r, r.readNegativeStrandFlag, tbe, r.cigarString, leftClipBases, rightClipBases, alignedBases
            )
        }
        return tbe
    }

    // we want to collect all these potential break ends and clean them up into distinct
    // locations, and then find supporting or nonsupporting evidence
    fun consolidatePotentialBreakEnds()
    {
        val combinedBreakEnds: MutableList<CandidateBreakEnd> = ArrayList()

        // first sort all the breakends
        // we want to sort it such that position is sorted last, and then we can use that
        // to compare if positions are close by, we merge them
        mCandidateBreakEnds.sortWith(
            Comparator.comparing(CandidateBreakEnd::key)
                .thenComparing(CandidateBreakEnd::chromosome)
                .thenComparingInt(CandidateBreakEnd::position)
        )

        // we want to merge all the ones that are the same
        var mergedBreakEnd: CandidateBreakEnd? = null
        for (tbe in mCandidateBreakEnds)
        {
            if (mergedBreakEnd != null && mergedBreakEnd.type === tbe.type &&
                mergedBreakEnd.chromosome == tbe.chromosome &&
                mergedBreakEnd.position == tbe.position
            )
            {
                mergedBreakEnd.splitReads.addAll(tbe.splitReads)
            }
            else
            {
                // make a new one
                mergedBreakEnd = tbe.copy()
                combinedBreakEnds.add(mergedBreakEnd)
            }
        }

        mCandidateBreakEnds.clear()

        // now we can go through this list and merge the ones that are close by to each other
        var maxSplitFragments = 0
        var currentBreakEnd: CandidateBreakEnd? = null
        for (tbe in combinedBreakEnds)
        {
            if (currentBreakEnd != null && currentBreakEnd.type === tbe.type &&
                currentBreakEnd.chromosome == tbe.chromosome &&
                abs(currentBreakEnd.position - tbe.position) <= BREAKEND_MERGE_DISTANCE
            )
            {
                if (tbe.splitReads.size > maxSplitFragments)
                {
                    // we can consolidate this one since it is near by. We count and see if number of
                    // supporting fragment is larger, and if it is then we use this position
                    currentBreakEnd.key.position = tbe.position
                    maxSplitFragments = tbe.splitReads.size
                }
                currentBreakEnd.splitReads.addAll(tbe.splitReads)
            } else
            {
                // make a new one
                currentBreakEnd = tbe.copy(splitReads = tbe.splitReads.toMutableList())
                maxSplitFragments = tbe.splitReads.size
                mCandidateBreakEnds.add(currentBreakEnd)
            }
        }
    }

    companion object
    {
        // where 1 end maps in the POLY-G region of LINC00486 (v38: chr2:32,916,190-32,916,630; v37: 2:33,141,260-33,141,700).
        fun isInExcludedBaseRegion(config: BreakEndParams, chromosome: String?, startPos: Int, endPos: Int): Boolean
        {
            var linc00486: ChrBaseRegion? = ExcludedRegions.getPolyGRegion(config.refGenomeVersion)

            if (linc00486 != null &&
                (linc00486.containsPosition(startPos) || linc00486.containsPosition(endPos)) &&
                ContigComparator.INSTANCE.compare(chromosome, linc00486.chromosome()) == 0)
            {
                return true
            }

            for (excludedRegion in config.excludedGenomeRegions)
            {
                if (startPos in excludedRegion.start() .. excludedRegion.end() &&
                    endPos in excludedRegion.start() .. excludedRegion.end() &&
                    ContigComparator.INSTANCE.compare(chromosome, excludedRegion.chromosome()) == 0)
                {
                    return true
                }
            }
            return false
        }
    }
}