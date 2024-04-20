package com.hartwig.hmftools.teal.breakend

import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeCoordinates
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion
import com.hartwig.hmftools.common.bam.CigarUtils
import com.hartwig.hmftools.common.bam.CigarUtils.*
import com.hartwig.hmftools.common.bam.SamRecordUtils
import com.hartwig.hmftools.teal.ReadGroup
import com.hartwig.hmftools.teal.TealUtils
import com.hartwig.hmftools.teal.util.TelomereMatcher
import htsjdk.samtools.SAMRecord
import org.apache.logging.log4j.LogManager
import kotlin.math.*

private const val MAX_DISTANCE_FROM_BREAKEND = 500

// split positions within 5 bases are count as the same
private const val SPLIT_POSITION_TOLERANCE = 5

// find all the disordant pairs given
// a list of candidate telomeric break ends
//
class BreakEndSupportCounter(refGenomeVersion: RefGenomeVersion, telomereMatchThreshold: Double)
{
    private val mLogger = LogManager.getLogger(javaClass)

    private val mTelomereMatchThreshold = telomereMatchThreshold
    private val mMinTelomereMatchLength = 12
    private val mRefGenomeCoordinates = when (refGenomeVersion) {
                                            RefGenomeVersion.V37 -> RefGenomeCoordinates.COORDS_37
                                            RefGenomeVersion.V38 -> RefGenomeCoordinates.COORDS_38 }

    private fun matchesTelomere(seq: String, gRich: Boolean): Boolean
    {
        return TelomereMatcher.matchesTelomere(seq, mTelomereMatchThreshold, mMinTelomereMatchLength, gRich)
    }

    // after we found all the candidate break ends, we want to process
    // the reads again and see if any discordant pairs support the split
    fun countBreakEndSupports(breakEnds: List<TelomericBreakEnd>, readGroups: Collection<ReadGroup>) : List<TelomericBreakEndSupport>
    {
        val breakEndSupports = breakEnds.map({ o -> TelomericBreakEndSupport(o) })

        for (rg in readGroups)
        {
            for (breakEnd: TelomericBreakEndSupport in breakEndSupports)
            {
                breakEndFragment(rg, breakEnd.key)?.apply({breakEnd.fragments.add(this)})
            }
        }

        // now we work out fragment length
        for (breakEnd: TelomericBreakEndSupport in breakEndSupports)
        {
            populateLongestTelomereSegment(breakEnd)
            populateLongestSplitReadAlign(breakEnd)
            populateDistanceFromTelomere(breakEnd)
        }

        return breakEndSupports
    }

    fun breakEndFragment(rg: ReadGroup, breakEnd: TelomericBreakEnd) : Fragment?
    {
        val alignedRead = findAlignedRead(rg, breakEnd)

        // we exclude duplicate reads
        if (alignedRead == null || alignedRead.duplicateReadFlag)
        {
            return null
        }

        val pairedRead = if (SamRecordUtils.firstInPair(alignedRead)) rg.secondOfPair else rg.firstOfPair

        if (pairedRead == null)
        {
            // shouldn't be here really, so we just do not classify this fragment
            return null
        }

        val alignedReadType = classifyAlignedRead(breakEnd, alignedRead)
        val pairedReadType = classifyPairedRead(breakEnd, pairedRead = pairedRead, alignedRead = alignedRead)

        if (alignedReadType != null)
        {
            mLogger.debug("breakend(${breakEnd}) readId(${rg.name}) aligned read type($alignedReadType) paired read type($pairedReadType)")
            // add it to the break end
            return Fragment(rg, alignedReadType = alignedReadType, pairedReadType = pairedReadType,
                alignedRead = alignedRead, pairedRead = pairedRead)
        }
        return null
    }

    // find the read that is aligned to the break point
    fun findAlignedRead(rg: ReadGroup, breakEnd: TelomericBreakEnd) : SAMRecord?
    {
        var alignedRead: SAMRecord? = null

        // first we find out which one is aligned close to the break point and which one
        // is the other read
        for (r in rg.allReads)
        {
            if (r.readUnmappedFlag)
            {
                continue
            }

            if (breakEnd.isRightTelomeric())
            {
                if (r.alignmentEnd <= breakEnd.position &&
                    r.alignmentEnd > breakEnd.position - MAX_DISTANCE_FROM_BREAKEND &&
                    r.referenceName == breakEnd.chromosome)
                {
                    // this is aligned read, but also compared to what we already found
                    // chose the one that is further to the right
                    if (alignedRead == null || r.alignmentEnd > alignedRead.alignmentEnd)
                    {
                        alignedRead = r
                    }
                }
            }
            else
            {
                // to support a new left telomere, we want to find a reverse reads to
                // be on the higher side of the break position, and the other
                // read to be fully telomeric
                if (r.alignmentStart >= breakEnd.position &&
                    r.alignmentStart < breakEnd.position + MAX_DISTANCE_FROM_BREAKEND &&
                    r.referenceName == breakEnd.chromosome)
                {
                    // this is aligned read, but also compared to what we already found
                    // chose the one that is further to the left
                    if (alignedRead == null || r.alignmentStart < alignedRead.alignmentStart)
                    {
                        alignedRead = r
                    }
                }
            }
        }

        // if the aligned read that we found is not facing the breakend and
        // does not cross it, then it is not interesting
        if (alignedRead != null)
        {
            if (breakEnd.isRightTelomeric() &&
                alignedRead.readNegativeStrandFlag &&
                alignedRead.alignmentEnd < breakEnd.position)
            {
                // for right telomere
                // reverse read faces away from breakend
                return null
            }

            if (!breakEnd.isRightTelomeric() &&
                !alignedRead.readNegativeStrandFlag &&
                alignedRead.alignmentStart > breakEnd.position)
            {
                // for left telomere
                // forward read faces away from breakend
                return null
            }
        }

        return alignedRead
    }

    fun classifyAlignedRead(breakEnd: TelomericBreakEnd, alignedRead: SAMRecord) : Fragment.AlignedReadType?
    {
        // check if aligned read is a split read
        var clipBases: String? = null

        if (breakEnd.isRightTelomeric() &&
            abs(alignedRead.alignmentEnd - breakEnd.position) <= SPLIT_POSITION_TOLERANCE)
        {
            clipBases = CigarUtils.rightSoftClipBases(alignedRead)
        }
        else if (!breakEnd.isRightTelomeric() &&
            abs(alignedRead.alignmentStart - breakEnd.position) <= SPLIT_POSITION_TOLERANCE)
        {
            clipBases = CigarUtils.leftSoftClipBases(alignedRead)
        }

        if (clipBases != null)
        {
            // now depending on whether the clipped bases are telomeric or not we assign types
            return if (matchesTelomere(clipBases, breakEnd.isGTelomeric()))
            {
                Fragment.AlignedReadType.SPLIT_READ_TELOMERIC
            }
            else
            {
                Fragment.AlignedReadType.SPLIT_READ_NOT_TELOMERIC
            }
        }

        // if this read spans across the break end, then we designate it as non split
        if ((alignedRead.alignmentStart + SPLIT_POSITION_TOLERANCE) < breakEnd.position &&
            (alignedRead.alignmentEnd - + SPLIT_POSITION_TOLERANCE) > breakEnd.position)
        {
            return Fragment.AlignedReadType.NOT_SPLIT_READ
        }

        // discordant pair. This usually means that the other
        if (breakEnd.isRightTelomeric() && alignedRead.alignmentEnd < breakEnd.position && !alignedRead.readNegativeStrandFlag ||
            !breakEnd.isRightTelomeric() && alignedRead.alignmentStart > breakEnd.position && alignedRead.readNegativeStrandFlag)
        {
            return Fragment.AlignedReadType.DISCORDANT_PAIR
        }

        return null
    }

    // We want to check the read group against the break end see if this supports
    // any of the following scenarios:
    // for right C telomere
    fun classifyPairedRead(breakEnd: TelomericBreakEnd,
                           pairedRead: SAMRecord,
                           alignedRead: SAMRecord) : Fragment.PairedReadType
    {
        // if the aligned read is not facing breakend then the pair read is not
        // interesting
        if (breakEnd.isRightTelomeric() == alignedRead.readNegativeStrandFlag)
        {
            return Fragment.PairedReadType.NOT_DISCORDANT
        }

        // now we need to make sure the other read sequence is represented properly
        // the negative strand flag is set based on where the paired read is aligned to in the ref
        // genome. However that alignment is not interesting for this purpose. We want to make sure the
        // paired read is represented such that it is pointing in the other direction of the aligned read.
        val pairReadString = if (alignedRead.readNegativeStrandFlag == pairedRead.readNegativeStrandFlag)
            TealUtils.reverseComplementSequence(pairedRead.readString) else pairedRead.readString

        // now we have filter down to the case where the aligned read is
        // aligned close to the correct side of the break point
        // we want to decide whether the read that is supposedly on the
        // other side of the break point is telomeric or not

        return if (matchesTelomere(pairReadString, breakEnd.isGTelomeric()))
                {
                    Fragment.PairedReadType.DISCORDANT_PAIR_TELOMERIC
                }
                else
                {
                    Fragment.PairedReadType.DISCORDANT_PAIR_NOT_TELOMERIC
                }
        // now we want to work out if this fragment contradict or support this break point
        /* TeloConfig.TE_LOGGER.info(
            "record({}) -v strand({}) breakend({}) cigar({}) leftClip({}) rightClip({}) aligned({})",
            r, r.readNegativeStrandFlag, tbe, r.cigarString, leftClipBases, rightClipBases, alignedBases
        )*/
    }

    private fun populateLongestTelomereSegment(breakEnd: TelomericBreakEndSupport)
    {
        // for all the reads, find the longest segment that matches our threshold
        val sequences = ArrayList<String>()

        for (frag in breakEnd.fragments)
        {
            // check the split reads
            val clipBases = if (breakEnd.type.isRightTelomeric())
                {
                    rightSoftClipBases(frag.alignedRead)
                }
                else
                {
                    leftSoftClipBases(frag.alignedRead)
                }
            if (clipBases != null)
            {
                sequences.add(clipBases)
            }

            if (frag.pairedReadType != Fragment.PairedReadType.NOT_DISCORDANT)
            {
                sequences.add(frag.pairedRead.readString)
            }
        }

        val telomereSegments = ArrayList<String>()

        // now we got all the sequences, we find the longest telomere match
        for (s in sequences)
        {
            val telomereMatch : TelomereMatcher.TelomereMatch? = if (breakEnd.type.isGTelomeric())
                {
                    TelomereMatcher.findGTelomereSegment(s, mTelomereMatchThreshold)
                }
                else
                {
                    TelomereMatcher.findCTelomereSegment(s, mTelomereMatchThreshold)
                }

            if (telomereMatch != null)
            {
                telomereSegments.add(telomereMatch.matchedSequence)
            }
        }

        // now find the longest one
        breakEnd.longestTelomereSegment = telomereSegments.maxByOrNull({ s -> s.length })
    }

    private fun populateLongestSplitReadAlign(breakEnd: TelomericBreakEndSupport)
    {
        // for all the reads, find the longest segment that matches our threshold
        for (frag in breakEnd.fragments)
        {
            // check the split reads
            if (frag.alignedReadType == Fragment.AlignedReadType.SPLIT_READ_TELOMERIC || frag.alignedReadType == Fragment.AlignedReadType.SPLIT_READ_NOT_TELOMERIC)
            {
                val splitRead = frag.alignedRead
                val alignedLength = splitRead.readLength - leftSoftClipLength(splitRead) - rightSoftClipLength(splitRead)
                breakEnd.longestSplitReadAlignLength = max(breakEnd.longestSplitReadAlignLength, alignedLength)
            }
        }
    }

    // we find the distance of the breakpoint to break end
    private fun populateDistanceFromTelomere(breakEnd: TelomericBreakEndSupport)
    {
        // find which ref genome to use
        var chromosomeLength: Int? = null

        try
        {
            chromosomeLength = mRefGenomeCoordinates.lengths()[HumanChromosome.fromString(breakEnd.chromosome)]
        }
        catch (_: IllegalArgumentException)
        {
        }

        if (chromosomeLength == null)
        {
            mLogger.debug("cannot find chromosome length for chromosome: {}", breakEnd.chromosome)
            return
        }

        val distance = min(breakEnd.position, chromosomeLength - breakEnd.position)

        breakEnd.distanceToTelomere = distance
    }
}
