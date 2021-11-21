package com.hartwig.hmftools.teal.breakend

import com.hartwig.hmftools.common.genome.chromosome.ContigComparator
import htsjdk.samtools.SAMRecord
import com.hartwig.hmftools.teal.TeloConfig
import com.hartwig.hmftools.common.samtools.SamRecordUtils
import com.hartwig.hmftools.teal.ReadGroup
import com.hartwig.hmftools.teal.TeloUtils
import com.hartwig.hmftools.teal.util.TelomereMatcher
import com.hartwig.hmftools.teal.TeloConstants
import java.util.Comparator
import kotlin.math.abs

// go through the bam file and find any potential
// split read sites
class SplitReadAnalyser(
    splitTelomereMatchThreshold: Double,
    alignedTelomereMatchThreshold: Double,
    fuzzyMatchDistance: Int,
    extraCandidateBreakEnds: Collection<TelomericBreakEndKey> = ArrayList())
{
    // we let some bases at the edge to be non telomeric
    private val mSplitTelomereMatchThreshold = splitTelomereMatchThreshold

    // I find using 0.8 can still get all the break ends gripss found.
    // using 0.7 and we start excluding some
    private val mAlignedTelomereMatchThreshold = alignedTelomereMatchThreshold
    private val mMinTelomereMatchLength = 12
    private val mFuzzyMatchDistance = fuzzyMatchDistance
    private var mBoundaryZone = 1

    private var mPotentialBreakEnds: MutableList<TelomericBreakEnd> = ArrayList()
    val potentialBreakEnds: List<TelomericBreakEnd> get() = mPotentialBreakEnds

    init
    {
        // we need to add the candidate breakends
        extraCandidateBreakEnds.forEach({ k -> mPotentialBreakEnds.add(TelomericBreakEnd(k)) })
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

        val tbe = findTelomericBreakEnd(r)
        if (tbe != null)
        {
            mPotentialBreakEnds.add(tbe)
            //TE_LOGGER.trace("record: {} breakend: {} cigar: {} readBases: {}", r, tbe, r.getCigarString(), r.getReadString());
        }
    }

    // look at soft clip site, and one part is mapped to a genome, and the
    // other part is telomeric
    // A right soft clip looks like 73S78M
    // A left soft clip looks like 78M32S
    // we want to work out whether one side is telomeric and the other side is not
    fun findTelomericBreakEnd(r: SAMRecord): TelomericBreakEnd?
    {
        // TE_LOGGER.trace("record: {}", r);
        if (isInExcludedBaseRegion(r.referenceName, r.alignmentStart, r.alignmentEnd))
        {
            return null
        }
        val readString = r.readString
        val leftClip = SamRecordUtils.leftSoftClip(r)
        val rightClip = SamRecordUtils.rightSoftClip(r)

        // now we work out if the clipped part is telomere
        val leftClipBases = readString.substring(0, leftClip)
        val rightClipBases = readString.substring(readString.length - rightClip)
        val alignedBases = readString.substring(leftClip, readString.length - rightClip)
        assert(leftClipBases.length == leftClip)
        assert(rightClipBases.length == rightClip)
        assert(leftClipBases.length + rightClipBases.length + alignedBases.length == readString.length)

        if (alignedBases.length <= mMinTelomereMatchLength)
        {
            // if the aligned section is too short to tell whether it is
            // even telomeric or not
            return null
        }

        if (TeloUtils.isPolyGC(alignedBases))
        {
            // aligned part must not be poly GC
            return null
        }

        TeloConfig.TE_LOGGER.trace("checking record: {}", r)

        var tbe: TelomericBreakEnd? = null
        if (leftClipBases.length >= mMinTelomereMatchLength + mBoundaryZone)
        {
            // we remove boundary zone
            val modLeftClipBases = leftClipBases.substring(0, leftClipBases.length - mBoundaryZone)
            if (matchesGTelomere(modLeftClipBases, mSplitTelomereMatchThreshold))
            {
                tbe = TelomericBreakEnd(TelomericBreakEndType.LEFT_G_TELOMERIC, r.referenceName, r.alignmentStart)
            }
            else if (matchesCTelomere(modLeftClipBases, mSplitTelomereMatchThreshold))
            {
                tbe = TelomericBreakEnd(TelomericBreakEndType.LEFT_C_TELOMERIC, r.referenceName, r.alignmentStart)
            }
        }
        if (tbe == null && rightClipBases.length >= mMinTelomereMatchLength + mBoundaryZone)
        {
            // we remove boundary zone
            val modRightClipBases = rightClipBases.substring(mBoundaryZone)
            if (matchesGTelomere(modRightClipBases, mSplitTelomereMatchThreshold))
            {
                tbe = TelomericBreakEnd(TelomericBreakEndType.RIGHT_G_TELOMERIC, r.referenceName, r.alignmentEnd)
            } else if (matchesCTelomere(modRightClipBases, mSplitTelomereMatchThreshold))
            {
                tbe = TelomericBreakEnd(TelomericBreakEndType.RIGHT_C_TELOMERIC, r.referenceName, r.alignmentEnd)
            }
        }

        if (tbe == null)
        {
            return null
        }

        // this could be potential telomeric region
        // next thing we need to check that the aligned region is not the same type of telomere
        when (tbe.type)
        {
            TelomericBreakEndType.LEFT_C_TELOMERIC, TelomericBreakEndType.RIGHT_C_TELOMERIC ->
            {
                if (matchesCTelomere(alignedBases, mAlignedTelomereMatchThreshold))
                {
                    // aligned bases is also  C telomeric, so this is no longer interesting
                    tbe = null
                }
            }
            TelomericBreakEndType.LEFT_G_TELOMERIC, TelomericBreakEndType.RIGHT_G_TELOMERIC ->
            {
                if (matchesGTelomere(alignedBases, mAlignedTelomereMatchThreshold))
                {
                    // aligned bases is also  C telomeric, so this is no longer interesting
                    tbe = null
                }
            }
        }

        if (tbe != null)
        {
            // add this to list of split reads
            tbe.telomericSplitReads.add(r)
            TeloConfig.TE_LOGGER.trace(
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
        val combinedBreakEnds: MutableList<TelomericBreakEnd> = ArrayList()

        // first sort all the breakends
        // we want to sort it such that position is sorted last, and then we can use that
        // to compare if positions are close by, we merge them
        mPotentialBreakEnds.sortWith(
            Comparator.comparing(TelomericBreakEnd::type)
                .thenComparing(TelomericBreakEnd::chromosome)
                .thenComparingInt(TelomericBreakEnd::position)
        )

        // first move we want to merge all the ones that are the same
        var currentBreakEnd: TelomericBreakEnd? = null
        for (tbe in potentialBreakEnds)
        {
            if (currentBreakEnd != null && currentBreakEnd.type === tbe.type &&
                currentBreakEnd.chromosome == tbe.chromosome &&
                currentBreakEnd.position == tbe.position)
            {
                // if they have same type and same chromosome, we check if they are nearby
                // we can combine them
                currentBreakEnd.telomericSplitReads.addAll(tbe.telomericSplitReads)
            }
            else
            {
                // make a new one
                currentBreakEnd = tbe.copy()
                combinedBreakEnds.add(currentBreakEnd)
            }
        }

        val consolidatedBreakEnds: MutableList<TelomericBreakEnd> = ArrayList()

        // now we can go through this list and merge the ones that are close by to each other
        var maxSplitFragments = 0
        currentBreakEnd = null
        for (tbe in combinedBreakEnds)
        {
            if (currentBreakEnd != null && currentBreakEnd.type === tbe.type &&
                currentBreakEnd.chromosome == tbe.chromosome &&
                abs(currentBreakEnd.position - tbe.position) <= mFuzzyMatchDistance
            )
            {
                if (tbe.telomericSplitReads.size > maxSplitFragments)
                {
                    // we can consolidate this one since it is near by. We count and see if number of
                    // supporting fragment is larger, and if it is then we use this position
                    currentBreakEnd.key.position = tbe.position
                    maxSplitFragments = tbe.telomericSplitReads.size
                }
                currentBreakEnd.telomericSplitReads.addAll(tbe.telomericSplitReads)
            } else
            {
                // make a new one
                currentBreakEnd = tbe.copy(telomericSplitReads = tbe.telomericSplitReads.toMutableList())
                maxSplitFragments = tbe.telomericSplitReads.size
                consolidatedBreakEnds.add(currentBreakEnd)
            }
        }

        // we just overwrite that list
        mPotentialBreakEnds = consolidatedBreakEnds
    }

    companion object
    {
        // where 1 end maps in the POLY-G region of LINC00486 (v38: chr2:32,916,190-32,916,630; v37: 2:33,141,260-33,141,700).
        fun isInExcludedBaseRegion(chromosome: String?, startPos: Int, endPos: Int): Boolean
        {
            for (excludedRegion in TeloConstants.EXCLUDED_BASE_REGIONS)
            {
                if (ContigComparator.INSTANCE.compare(chromosome, excludedRegion.chromosome()) == 0 &&
                    (excludedRegion.containsPosition(startPos) || excludedRegion.containsPosition(endPos))
                )
                {
                    return true
                }
            }
            return false
        }
    }
}