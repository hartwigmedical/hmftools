package com.hartwig.hmftools.esvee.prep;

import static java.lang.Math.abs;
import static java.lang.Math.min;

import static com.hartwig.hmftools.common.bam.SamRecordUtils.ALIGNMENT_SCORE_ATTRIBUTE;
import static com.hartwig.hmftools.common.bam.SamRecordUtils.XS_ATTRIBUTE;
import static com.hartwig.hmftools.common.region.BaseRegion.positionWithin;
import static com.hartwig.hmftools.esvee.assembly.AssemblyConfig.SV_LOGGER;
import static com.hartwig.hmftools.esvee.common.CommonUtils.belowMinQual;
import static com.hartwig.hmftools.esvee.common.SvConstants.MAX_JUNC_POSITION_DIFF;
import static com.hartwig.hmftools.esvee.prep.PrepConstants.MAX_HIGH_QUAL_BASE_MISMATCHES;
import static com.hartwig.hmftools.esvee.prep.PrepConstants.MIN_ALIGNMENT_SCORE_DIFF;
import static com.hartwig.hmftools.esvee.prep.PrepConstants.MIN_CALC_ALIGNMENT_LOWER_SCORE;
import static com.hartwig.hmftools.esvee.prep.PrepConstants.MIN_CALC_ALIGNMENT_SCORE;
import static com.hartwig.hmftools.esvee.prep.PrepConstants.MIN_EXACT_BASE_PERC;
import static com.hartwig.hmftools.esvee.prep.ReadFilters.calcRepeatTrimmedAlignmentScore;
import static com.hartwig.hmftools.esvee.prep.ReadFilters.isChimericRead;
import static com.hartwig.hmftools.esvee.prep.types.FragmentData.unclippedPosition;
import static com.hartwig.hmftools.esvee.prep.types.ReadGroupStatus.DUPLICATE;

import java.util.Collection;
import java.util.List;
import java.util.Map;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.genome.region.Orientation;
import com.hartwig.hmftools.esvee.common.ReadIdTrimmer;
import com.hartwig.hmftools.esvee.prep.types.FragmentData;
import com.hartwig.hmftools.esvee.prep.types.JunctionData;
import com.hartwig.hmftools.esvee.prep.types.PrepRead;
import com.hartwig.hmftools.esvee.prep.types.ReadFilterConfig;
import com.hartwig.hmftools.esvee.prep.types.ReadGroup;
import com.hartwig.hmftools.esvee.prep.types.ReadType;

public final class JunctionUtils
{
    protected static boolean readWithinJunctionRange(final PrepRead read, final JunctionData junctionData, int distanceBuffer)
    {
        if(abs(read.AlignmentEnd - junctionData.Position) <= distanceBuffer)
            return true;

        if(abs(read.AlignmentStart - junctionData.Position) <= distanceBuffer)
            return true;

        return false;
    }

    protected static final int INVALID_JUNC_INDEX = -1;
    protected static final int SIMPLE_SEARCH_COUNT = 10;

    @VisibleForTesting
    protected static int findJunctionIndex(final List<JunctionData> junctions, int position, int simpleSearchCount)
    {
        // returns the index of a matched junction position, otherwise the preceding index, otherwise invalid if this is lower
        if(junctions.isEmpty() || position < junctions.get(0).Position)
            return INVALID_JUNC_INDEX;

        int junctionCount = junctions.size();
        int maxIndex = junctionCount - 1;

        if(position >= junctions.get(junctionCount - 1).Position)
            return maxIndex;

        if(junctionCount <= simpleSearchCount)
        {
            for(int index = 0; index < junctionCount; ++index)
            {
                JunctionData junctionData = junctions.get(index);

                if(junctionData.Position < position)
                    continue;

                if(junctionData.Position > position)
                    return index - 1;

                return index;
            }

            return maxIndex;
        }

        // binary search on junctions for a larger collection
        int currentIndex = junctions.size() / 2;
        int lowerIndex = 0;
        int upperIndex = maxIndex;

        int iterations = 0;

        while(true)
        {
            JunctionData junctionData = junctions.get(currentIndex);

            if(junctionData.Position == position)
                return currentIndex;

            if(upperIndex == lowerIndex + 1)
                break;

            if(position < junctionData.Position)
            {
                // need to search lower
                if(currentIndex == 0)
                    break;

                if(currentIndex == 1)
                {
                    if(lowerIndex == currentIndex)
                        break;

                    upperIndex = 1;
                    currentIndex = 0;
                }
                else
                {
                    upperIndex = currentIndex;
                    currentIndex = (lowerIndex + upperIndex) / 2;
                }
            }
            else if(position > junctionData.Position)
            {
                // search higher
                if(currentIndex == maxIndex)
                    break;

                if(currentIndex == maxIndex - 1)
                {
                    if(upperIndex == currentIndex)
                        break;

                    lowerIndex = currentIndex;
                    currentIndex = maxIndex;
                }
                else
                {
                    lowerIndex = currentIndex;
                    currentIndex = (lowerIndex + upperIndex) / 2;
                }
            }
            else
            {
                break;
            }

            ++iterations;

            if(iterations > 100)
            {
                SV_LOGGER.warn("junction index search iterations({}}) junctions({}) index(cur={} low={} high={})",
                        iterations, junctions.size(), currentIndex, lowerIndex, upperIndex);
                break;
            }
        }

        JunctionData junctionData = junctions.get(currentIndex);

        // if not found, return the index of the junction immediately below this range
        if(position > junctionData.Position)
            return currentIndex;
        else
            return currentIndex - 1;
    }

    protected static int findJunctionMatchIndex(final List<JunctionData> junctions, int position, final Orientation orientation, int juncIndex)
    {
        // finds an exact match on position and orientation around the current index
        for(int i = -1; i <= 1; ++i)
        {
            int index = juncIndex + i;

            if(index >= 0 && index < junctions.size())
            {
                JunctionData junctionData = junctions.get(index);

                if(junctionData.Position == position && junctionData.Orient == orientation)
                    return index;
            }
        }

        return INVALID_JUNC_INDEX;
    }

    public static boolean hasOtherJunctionSupport(
            final PrepRead read, final JunctionData junctionData, final ReadFilterConfig filterConfig)
    {
        int unclippedStart = read.UnclippedStart;
        int unclippedEnd = read.UnclippedEnd;

        // first check for a read crossing the junction
        if(positionWithin(junctionData.Position, unclippedStart, unclippedEnd))
        {
            // correct side of the junction
            int junctionDistance = 0;

            if(junctionData.isForward())
            {
                junctionDistance = min(abs(unclippedEnd - junctionData.Position), abs(read.AlignmentEnd - junctionData.Position));
            }
            else
            {
                junctionDistance = min(abs(unclippedStart - junctionData.Position), abs(read.AlignmentStart - junctionData.Position));
            }

            // any soft-clipping on the correct side if close to the junction
            if(junctionDistance <= MAX_JUNC_POSITION_DIFF)
            {
                if(junctionData.isForward() && read.isRightClipped())
                    return true;

                if(junctionData.isReverse() && read.isLeftClipped())
                    return true;
            }

            return false;
        }

        // otherwise can be distant if discordant and with an orientation cross the junction
        int junctionDistance = 0;

        if(junctionData.Orient != read.orientation())
            return false;

        if(junctionData.isForward())
        {
            if(read.AlignmentEnd > junctionData.Position)
                return false;

            junctionDistance = abs(read.AlignmentEnd - junctionData.Position);
        }
        else
        {
            if(read.AlignmentStart < junctionData.Position) //  || abs(read.AlignmentEnd - junctionData.Position) > filterConfig.maxSupportingFragmentDistance()
                return false;

            junctionDistance = abs(read.AlignmentStart - junctionData.Position);
        }

        if(junctionDistance <= filterConfig.maxSupportingFragmentDistance())
            return isChimericRead(read.record(), filterConfig);

        return false;
    }

    public static boolean hasExactJunctionSupport(
            final PrepRead read, final JunctionData junctionData, final ReadFilterConfig filterConfig)
    {
        // for a read to be classified as exact support it needs to meet the following criteria:
        // a) soft or hard-clipped at exactly the same base as the junction, otherwise
        // b) soft-clipped before or after the junction with:
        // - the read's ref/SC bases matching any overlapping junction ref/SC bases
        // - allowing for 1 high-qual mismatch
        // - ignoring low-qual mismatches
        // - requiring > 25% of all bases to match regardless of qual

        int readStart = read.AlignmentStart;
        int readEnd = read.AlignmentEnd;
        int readLength = read.readBases().length();

        PrepRead juncRead = junctionData.topJunctionRead();

        if(juncRead == null)
            return false;

        int juncReadLength = juncRead.readBases().length();
        int juncReadCigarCount = juncRead.cigar().getCigarElements().size();

        if(juncReadCigarCount < 2)
            return false;

        // the comparison of bases between this read and the designed junction read will start at the inner of the 2 read's alignments
        int readIndexStart, juncReadIndexStart;

        if(junctionData.isForward())
        {
            if(!read.isRightClipped())
                return false;

            if(readEnd == junctionData.Position) // an exact match is checked no further in prep
                return true;

            if(abs(readEnd - junctionData.Position) > MAX_JUNC_POSITION_DIFF) // aligned too far from the junction
                return false;

            // read and soft-clipped position must straddle the junction
            if(!(readStart < junctionData.Position && read.UnclippedEnd > junctionData.Position))
                return false;

            int juncPosDiff = juncRead.AlignmentEnd - read.AlignmentEnd;

            if(juncPosDiff < 0)
            {
                // junction read ends earlier - comparison will start from it's first soft-clip base and include some of the read's aligned bases
                readIndexStart = readLength - read.rightClipLength() + juncPosDiff; // note taking off a -ve adjustment
                juncReadIndexStart = juncReadLength - juncRead.rightClipLength();
            }
            else
            {
                // comparison will include some of the junction read's ref bases
                readIndexStart = readLength - read.rightClipLength();
                juncReadIndexStart = juncReadLength - juncRead.rightClipLength() - juncPosDiff;
            }
        }
        else
        {
            if(!read.isLeftClipped())
                return false;

            if(readStart == junctionData.Position)
                return true;

            if(abs(readStart - junctionData.Position) > MAX_JUNC_POSITION_DIFF)
                return false;

            if(!(readEnd > junctionData.Position && read.UnclippedStart < junctionData.Position))
                return false;

            int juncPosDiff = juncRead.AlignmentStart - read.AlignmentStart;

            if(juncPosDiff < 0)
            {
                // junction read starts earlier
                readIndexStart = read.leftClipLength() - 1;
                juncReadIndexStart = juncRead.leftClipLength() - 1 + (-juncPosDiff); // adding a -ve adjustment
            }
            else
            {
                readIndexStart = read.leftClipLength() - 1 + juncPosDiff;
                juncReadIndexStart = juncRead.leftClipLength() - 1;
            }
        }

        // work out bases to compare with the designated junction read
        boolean searchDown = junctionData.isReverse();
        return basesMeetExactSupport(juncRead, read, juncReadIndexStart, readIndexStart, searchDown);
    }

    private static boolean basesMeetExactSupport(
            final PrepRead juncRead, final PrepRead read, int juncReadIndex, int readIndex, boolean searchDown)
    {
        int highQualMismatches = 0;
        int baseMatches = 0;

        int readIndexMax = read.readBases().length() - 1;
        int juncReadIndexMax = juncRead.readBases().length() - 1;
        int basesCompared = 0;

        while(juncReadIndex >= 0 && readIndex >= 0 && juncReadIndex <= juncReadIndexMax && readIndex <= readIndexMax)
        {
            char readBase = read.readBases().charAt(readIndex);
            char juncReadBase = juncRead.readBases().charAt(juncReadIndex);
            ++basesCompared;

            if(readBase == juncReadBase)
            {
                ++baseMatches;
            }
            else
            {
                if(!belowMinQual(read.baseQualities()[readIndex]) && !belowMinQual(juncRead.baseQualities()[juncReadIndex]))
                {
                    ++highQualMismatches;

                    if(highQualMismatches > MAX_HIGH_QUAL_BASE_MISMATCHES)
                        return false;
                }
            }

            juncReadIndex += searchDown ? -1 : 1;
            readIndex += searchDown ? -1 : 1;
        }

        double baseMatchPerc = baseMatches / (double)basesCompared;
        return baseMatchPerc >= MIN_EXACT_BASE_PERC;
    }

    public static boolean hasWellAnchoredRead(final JunctionData junctionData, final ReadFilterConfig filterConfig)
    {
        Collection<PrepRead> junctionReads = junctionData.readTypeReads().get(ReadType.JUNCTION);

        if(junctionReads.stream().anyMatch(x -> x.hasLineTail()))
            return true;

        for(PrepRead read : junctionReads)
        {
            double adjustedAlignScore = calcRepeatTrimmedAlignmentScore(read, MIN_CALC_ALIGNMENT_SCORE, true);

            if(adjustedAlignScore < 0)
                continue;

            if(adjustedAlignScore >= MIN_CALC_ALIGNMENT_SCORE)
                return true;

            if(adjustedAlignScore >= MIN_CALC_ALIGNMENT_LOWER_SCORE && aboveAlignedScoreDifference(read, filterConfig.MinAlignmentBases))
                return true;
        }

        return false;
    }

    private static boolean aboveAlignedScoreDifference(final PrepRead read, int minAlignScore)
    {
        Integer asScore = read.record().getIntegerAttribute(ALIGNMENT_SCORE_ATTRIBUTE);

        if(asScore == null || asScore < minAlignScore)
            return false;

        Integer xsScore = read.record().getIntegerAttribute(XS_ATTRIBUTE);

        return xsScore != null && asScore - xsScore >= MIN_ALIGNMENT_SCORE_DIFF;
    }

    public static int markSupplementaryDuplicates(
            final Map<String,ReadGroup> readGroupMap, final ReadIdTrimmer readIdTrimmer, int permittedPositionDiff)
    {
        Map<Integer,List<PrepRead>> initialPositionMap = Maps.newHashMap();

        for(ReadGroup readGroup : readGroupMap.values())
        {
            if(readGroup.reads().stream().noneMatch(x -> x.hasSuppAlignment()))
                continue;

            if(readGroup.reads().stream().noneMatch(x -> x.readType() == ReadType.JUNCTION || x.readType() == ReadType.CANDIDATE_SUPPORT))
                continue;

            // require paired reads with a supplementary read
            for(PrepRead read : readGroup.reads())
            {
                if(!read.hasSuppAlignment())
                    continue;

                if(read.isPaired() && !read.isMateMapped())
                    continue;

                if(read.readType() != ReadType.JUNCTION && read.readType() != ReadType.CANDIDATE_SUPPORT)
                    continue;

                int unclippedPosition = unclippedPosition(read);

                List<PrepRead> matchingGroups = null;

                if(permittedPositionDiff == 0)
                {
                    matchingGroups = initialPositionMap.get(unclippedPosition);
                }
                else
                {
                    for(Map.Entry<Integer,List<PrepRead>> entry : initialPositionMap.entrySet())
                    {
                        if(abs(entry.getKey() - unclippedPosition) <= permittedPositionDiff)
                        {
                            matchingGroups = entry.getValue();
                            break;
                        }
                    }
                }

                if(matchingGroups == null)
                {
                    matchingGroups = Lists.newArrayList();
                    initialPositionMap.put(unclippedPosition, matchingGroups);
                }

                matchingGroups.add(read);
                break;
            }
        }

        int duplicateFragmentCount = 0;

        for(List<PrepRead> reads : initialPositionMap.values())
        {
            if(reads.size() < 2)
                continue;

            boolean hasSupp = false;
            boolean hasPrimary = false;

            List<FragmentData> fragments = Lists.newArrayListWithCapacity(reads.size());

            for(PrepRead read : reads)
            {
                fragments.add(new FragmentData(read));

                if(read.isSupplementaryAlignment())
                    hasSupp = true;
                else
                    hasPrimary = true;
            }

            if(!hasSupp || !hasPrimary)
                continue;

            // mark any group where the supplementary is a duplicate and the lower (in coord terms) of the two spanning reads
            for(int i = 0; i < reads.size() - 1; ++i)
            {
                FragmentData firstFrag = fragments.get(i);

                for(int j = i + 1; j < reads.size(); ++j)
                {
                    FragmentData nextFrag = fragments.get(j);

                    if(firstFrag.IsPrimary == nextFrag.IsPrimary)
                        continue;

                    if(permittedPositionDiff == 0)
                    {
                        if(!firstFrag.matches(nextFrag))
                            continue;
                    }
                    else
                    {
                        if(!firstFrag.withinPositionRange(nextFrag, permittedPositionDiff))
                            continue;
                    }

                    // the supp will be removed if it is the lower coordinate of the 2, or its supplementary is likewise
                    // but favour consensus reads over non-consensus
                    boolean markNextAsDuplicate;

                    boolean firstIsConsensus = firstFrag.isConsensus();

                    if(firstIsConsensus != nextFrag.isConsensus())
                    {
                        markNextAsDuplicate = firstIsConsensus;
                    }
                    else
                    {
                        markNextAsDuplicate = (firstFrag.IsPrimary == firstFrag.ReadIsLowerVsSuppData);
                    }

                    String duplicateReadId = readIdTrimmer.trim(markNextAsDuplicate ? nextFrag.Read.id() : firstFrag.Read.id());

                    SV_LOGGER.trace("duplicate fragment({}) vs other({})",
                            markNextAsDuplicate ? nextFrag : firstFrag, markNextAsDuplicate ? firstFrag : nextFrag);

                    ReadGroup duplicateGroup = readGroupMap.get(duplicateReadId);

                    if(duplicateGroup == null)
                    {
                        SV_LOGGER.error("duplicateReadId({}) group not found from read({}), trimmerEnabled({})",
                                duplicateReadId, markNextAsDuplicate ? nextFrag.Read.id() : firstFrag.Read.id(), readIdTrimmer.enabled());
                        return duplicateFragmentCount;
                    }

                    duplicateGroup.setGroupState(DUPLICATE);
                    ++duplicateFragmentCount;
                }
            }
        }

        return duplicateFragmentCount;
    }
}
