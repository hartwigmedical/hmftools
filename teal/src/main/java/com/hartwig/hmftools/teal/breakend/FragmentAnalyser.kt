package com.hartwig.hmftools.teal.breakend

import com.hartwig.hmftools.common.samtools.SamRecordUtils
import com.hartwig.hmftools.teal.ReadGroup
import com.hartwig.hmftools.teal.TeloUtils
import com.hartwig.hmftools.teal.util.TelomereMatcher
import htsjdk.samtools.SAMRecord
import org.apache.logging.log4j.LogManager
import kotlin.math.abs

private const val MAX_DISTANCE_FROM_BREAKEND = 500

// find all the disordant pairs given
// a list of candidate telomeric break ends
//
class FragmentAnalyser(breakEnds: List<TelomericBreakEnd>, fuzzyMatchDistance: Int)
{
    private val LOGGER = LogManager.getLogger(javaClass)

    private data class ReadAlignment(val alignedReads: List<SAMRecord> = ArrayList(), val otherReads: List<SAMRecord> = ArrayList())

    // TODO: this should be linked to fragment size
    private val mTelomereMatchThreshold = 0.9
    private val mMinTelomereMatchLength = 12
    private val mBreakEnds = breakEnds
    private val mFuzzyMatchDistance = fuzzyMatchDistance

    val breakEnds: List<TelomericBreakEnd> get() = mBreakEnds

    private fun matchesGTelomere(seq: String): Boolean
    {
        return TelomereMatcher.matchesGTelomere(seq, mTelomereMatchThreshold, mMinTelomereMatchLength)
    }

    private fun matchesCTelomere(seq: String): Boolean
    {
        return TelomereMatcher.matchesCTelomere(seq, mTelomereMatchThreshold, mMinTelomereMatchLength)
    }

    // after we found all the candidate break ends, we want to process
    // the reads again and see if any discordant pairs support the split
    fun processReadGroups(readGroups: Collection<ReadGroup>)
    {
        for (rg in readGroups)
        {
            for (breakEnd: TelomericBreakEnd in mBreakEnds)
            {
                checkAgainstBreakEnd(rg, breakEnd)
            }
        }
    }

    fun checkAgainstBreakEnd(rg: ReadGroup, breakEnd: TelomericBreakEnd)
    {
        val facingBreakReads : Collection<SAMRecord> = findReadFacingBreak(rg, breakEnd)

        val splitReads = findSplitReads(rg, breakEnd)

        if (facingBreakReads.isEmpty() && splitReads.isEmpty())
        {
            return
        }

        var fragType : Fragment.Type? = if (facingBreakReads.isNotEmpty()) checkForDiscordantPair(rg, breakEnd, facingBreakReads) else null

        if (splitReads.isNotEmpty())
        {
            when (fragType)
            {
                null, Fragment.Type.DISCORDANT_PAIR -> fragType = Fragment.Type.SPLIT_READ
                Fragment.Type.CONTRA_DISCORDANT_PAIR -> fragType = Fragment.Type.CONTRA_SPLIT_READ
                else -> {}
            }
        }

        if (fragType != null)
        {
            // add it to the break end
            breakEnd.fragments.add(Fragment(fragType, rg, splitReads = splitReads, facingBreakReads = facingBreakReads))
        }
    }

    // first stage we want to actually classify the reads according to this breakend
    private fun findReadFacingBreak(rg: ReadGroup, breakEnd: TelomericBreakEnd) : Collection<SAMRecord>
    {
        val facingBreakReads = ArrayList<SAMRecord>()

        val isRightTelomeric = breakEnd.type.isRightTelomeric()

        // first we find out which one is aligned close to the break point and which one
        // is the other read
        for (r in rg.allReads)
        {
            if (isRightTelomeric)
            {
                // to support a new right telomere, we want to find a forward read to
                // map to left side of the break position, and the other
                // read to be fully telomeric
                if (!r.readUnmappedFlag &&
                    !r.readNegativeStrandFlag &&
                    r.alignmentEnd <= breakEnd.position &&
                    r.alignmentEnd > breakEnd.position - MAX_DISTANCE_FROM_BREAKEND &&
                    r.referenceName == breakEnd.chromosome)
                {
                    facingBreakReads.add(r)
                }
            }
            else
            {
                // to support a new left telomere, we want to find a reverse reads to
                // be on the higher side of the break position, and the other
                // read to be fully telomeric
                if (!r.readUnmappedFlag &&
                    r.readNegativeStrandFlag &&
                    r.alignmentStart >= breakEnd.position &&
                    r.alignmentStart < breakEnd.position + MAX_DISTANCE_FROM_BREAKEND &&
                    r.referenceName == breakEnd.chromosome)
                {
                    facingBreakReads.add(r)
                }
            }
        }

        if (facingBreakReads.isEmpty())
        {
            return emptyList()
        }

        // special case, if we got more than 1 aligned reads then we need to work out if the reads are actually
        // overlapping. This would be a sign that there is a repeated section inserted.
        // since we know both reads are aligned to same location we can rely on first of pair flag
        if (facingBreakReads.size == 2)
        {
            if (facingBreakReads[0].firstOfPairFlag == facingBreakReads[1].firstOfPairFlag)
            {
                // we got same read that have multiple alignment to the same location. Just choose the one
                // with the better score???
                LOGGER.debug("breakend({}) multiple alignment of same read: {} {}", breakEnd.key, facingBreakReads[0], facingBreakReads[1])
            }
            else
            {
                // this is very unusual, we got both reads pointing at the same direction aligned at the same location
                // log it for now
                LOGGER.debug("breakend({}) potential short inversion reads: {} {}", breakEnd.key, facingBreakReads[0], facingBreakReads[1])
            }
        }

        // if we only have one aligned read then we need to work out whether the other read is on the other
        // side of the breakend

        return facingBreakReads
    }

    // We want to check the read group against the break end see if this supports
    // any of the following scenarios:
    // for right C telomere
    fun checkForDiscordantPair(rg: ReadGroup, breakEnd: TelomericBreakEnd, readsFacingBreak: Collection<SAMRecord>) : Fragment.Type?
    {
        if (readsFacingBreak.isEmpty())
        {
            // shouldn't be here
            return null
        }

        if (readsFacingBreak.size > 1)
        {
            // if there are multiple reads aligned to same location, we want to see if both reads are aligned, if so we
            // cannot count this as discordant pair
            if (readsFacingBreak.any({ r -> r.firstOfPairFlag }) && readsFacingBreak.any({ r -> !r.firstOfPairFlag }))
            {
                // both reads are aligned to the break point
                return null
            }
        }

        val alignedRead: SAMRecord = readsFacingBreak.first()
        val otherRead = if (alignedRead == rg.firstOfPair) rg.secondOfPair else rg.firstOfPair

        // now we need to make sure the other read sequence is represented properly
        // the reason is that we might not be aligning them correctly, so we want to make
        // sure that read reversed flag is ok
        val otherReadString = if (alignedRead.readNegativeStrandFlag == otherRead.readNegativeStrandFlag)
            TeloUtils.reverseComplementSequence(otherRead.readString) else otherRead.readString

        // now we have filter down to the case where the aligned read is
        // aligned close to the correct side of the break point
        // we want to decide whether the read that is supposedly on the
        // other side of the break point is telomeric or not
        when (breakEnd.type)
        {
            TelomericBreakEndType.RIGHT_G_TELOMERIC, TelomericBreakEndType.LEFT_G_TELOMERIC ->
            {
                if (matchesGTelomere(otherReadString))
                {
                    // supporting
                    return Fragment.Type.DISCORDANT_PAIR
                }
                else
                {   // contradicting
                    return Fragment.Type.CONTRA_DISCORDANT_PAIR
                }
            }
            TelomericBreakEndType.RIGHT_C_TELOMERIC, TelomericBreakEndType.LEFT_C_TELOMERIC ->
            {
                if (matchesCTelomere(otherReadString))
                {
                    // supporting
                    return Fragment.Type.DISCORDANT_PAIR
                }
                else
                {   // contradicting
                    return Fragment.Type.CONTRA_DISCORDANT_PAIR
                }
            }
        }
        // now we want to work out if this fragment contradict or support this break point
        /* TeloConfig.TE_LOGGER.info(
            "record({}) -v strand({}) breakend({}) cigar({}) leftClip({}) rightClip({}) aligned({})",
            r, r.readNegativeStrandFlag, tbe, r.cigarString, leftClipBases, rightClipBases, alignedBases
        )*/
    }

    // We want to check the read group against the break end see if this supports
    // any of the following scenarios:
    // for right C telomere
    private fun findSplitReads(rg: ReadGroup, breakEnd: TelomericBreakEnd) : Collection<SAMRecord>
    {
        val splitReads = ArrayList<SAMRecord>()

        // check if any of the reads are split on this break end
        for (r in rg.allReads)
        {
            if (r.readUnmappedFlag)
            {
                continue
            }

            if (breakEnd.isRightTelomeric() &&
                abs(r.alignmentEnd - breakEnd.position) < mFuzzyMatchDistance &&
                r.referenceName == breakEnd.chromosome)
            {
                if (SamRecordUtils.rightSoftClip(r) != 0)
                {
                    splitReads.add(r)
                }
            }
            else if (!breakEnd.isRightTelomeric() &&
                abs(r.alignmentStart - breakEnd.position) < mFuzzyMatchDistance &&
                r.referenceName == breakEnd.chromosome)
            {
                if (SamRecordUtils.leftSoftClip(r) != 0)
                {
                    splitReads.add(r)
                }
            }
        }

        return splitReads
    }

}
