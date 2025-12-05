package com.hartwig.hmftools.esvee.prep;

import static java.lang.Math.abs;
import static java.lang.Math.max;
import static java.lang.Math.min;

import static com.hartwig.hmftools.common.bam.SamRecordUtils.ALIGNMENT_SCORE_ATTRIBUTE;
import static com.hartwig.hmftools.common.bam.SamRecordUtils.XS_ATTRIBUTE;
import static com.hartwig.hmftools.common.region.BaseRegion.positionWithin;
import static com.hartwig.hmftools.esvee.assembly.AssemblyConfig.SV_LOGGER;
import static com.hartwig.hmftools.esvee.common.CommonUtils.belowMinQual;
import static com.hartwig.hmftools.esvee.prep.PrepConstants.MAX_HIGH_QUAL_BASE_MISMATCHES;
import static com.hartwig.hmftools.esvee.prep.PrepConstants.MIN_ALIGNMENT_SCORE_DIFF;
import static com.hartwig.hmftools.esvee.prep.PrepConstants.MIN_CALC_ALIGNMENT_LOWER_SCORE;
import static com.hartwig.hmftools.esvee.prep.PrepConstants.MIN_CALC_ALIGNMENT_SCORE;
import static com.hartwig.hmftools.esvee.prep.PrepConstants.MIN_EXACT_BASE_PERC;
import static com.hartwig.hmftools.esvee.prep.ReadFilters.calcRepeatTrimmedAlignmentScore;
import static com.hartwig.hmftools.esvee.prep.ReadFilters.isChimericRead;
import static com.hartwig.hmftools.esvee.prep.types.FragmentData.unclippedPosition;
import static com.hartwig.hmftools.esvee.prep.types.ReadGroupStatus.DUPLICATE;

import static htsjdk.samtools.CigarOperator.M;
import static htsjdk.samtools.CigarOperator.S;

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

import htsjdk.samtools.CigarElement;

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
            if(junctionDistance <= filterConfig.MinSupportingReadDistance)
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
        boolean leftSoftClipped = read.cigar().isLeftClipped();
        boolean rightSoftClipped = read.cigar().isRightClipped();

        if(!leftSoftClipped && !rightSoftClipped)
            return false;

        // for a read to be classified as exact support it needs to meet the following criteria:
        // a) soft or hard-clipped at exactly the same base as the junction
        // b) soft-clipped before or after the junction with:
        // - the read's ref/SC bases matching any overlapping junction ref/SC bases
        // - allowing for 1 high-qual mismatch
        // - ignoring low-qual mismatches
        // - requiring > 25% of all bases to match

        final PrepRead juncRead = junctionData.topJunctionRead();

        int readLength = read.readBases().length();

        if(junctionData.isForward())
        {
            if(!rightSoftClipped)
                return false;

            int readRightPos = read.AlignmentEnd;

            if(readRightPos == junctionData.Position)
                return true;

            if(juncRead == null)
                return false;

            // within 50 bases with exact sequence match in between the soft clip locations
            if(abs(readRightPos - junctionData.Position) > filterConfig.MinSupportingReadDistance)
                return false;

            int scLength = 0;
            int firstMatchLength = 0;

            for(int i = read.cigar().getCigarElements().size() - 1 ; i >= 0; --i)
            {
                CigarElement element = read.cigar().getCigarElements().get(i);

                if(element.getOperator() == S)
                {
                    scLength = element.getLength();
                }
                else if(element.getOperator() == M)
                {
                    firstMatchLength = element.getLength();
                    break;
                }
            }

            // must also overlap the junction
            if(read.AlignmentStart > junctionData.Position || readRightPos + scLength < junctionData.Position)
                return false;

            int readEndPosIndex = readLength - scLength - 1;

            int juncReadLength = juncRead.readBases().length();
            int juncReadScLength = juncRead.cigar().getLastCigarElement().getLength();
            int juncReadEndPosIndex = juncReadLength - juncReadScLength - 1;
            int endPosDiff = juncRead.AlignmentEnd - readRightPos;

            int junctionReadOffset = juncReadEndPosIndex - readEndPosIndex - endPosDiff;

            // test all overlapping bases - either from ref or soft-clip bases
            int startIndex = readLength - scLength - min(max(read.AlignmentEnd - junctionData.Position, 0), firstMatchLength);

            if(startIndex < 0)
                return false;

            int highQualMismatches = 0;
            int baseMatches = 0;
            for(int i = startIndex; i < readLength; ++i)
            {
                char readBase = read.readBases().charAt(i);

                int juncIndex = i + junctionReadOffset;
                if(juncIndex < 0 || juncIndex >= juncReadLength)
                    return false;

                char juncReadBase = juncRead.readBases().charAt(juncIndex);

                if(readBase == juncReadBase)
                {
                    ++baseMatches;
                    continue;
                }

                if(belowMinQual(read.baseQualities()[i]) || belowMinQual(juncRead.baseQualities()[juncIndex]))
                    continue;

                ++highQualMismatches;

                if(highQualMismatches > MAX_HIGH_QUAL_BASE_MISMATCHES)
                    return false;
            }

            double baseMatchPerc = baseMatches / (double)(readLength - startIndex);
            return baseMatchPerc > MIN_EXACT_BASE_PERC;
        }
        else
        {
            // negative orientation
            if(!leftSoftClipped)
                return false;

            int readLeftPos = read.AlignmentStart;

            if(readLeftPos == junctionData.Position)
                return true;

            if(juncRead == null)
                return false;

            // within 50 bases with exact sequence match in between the soft clip locations
            if(abs(readLeftPos - junctionData.Position) > filterConfig.MinSupportingReadDistance)
                return false;

            // test for a base match for the read's soft-clipped bases, allow for low-qual matches

            // read: SC length -> start position
            // junc: SC length -> start position
            // junc read index = sc length diff - position diff

            int scLength = 0;
            int firstMatchLength = 0;

            for(CigarElement element : read.cigar().getCigarElements())
            {
                if(element.getOperator() == S)
                {
                    scLength = element.getLength();
                }
                else if(element.getOperator() == M)
                {
                    firstMatchLength = element.getLength();
                    break;
                }
            }

            if(read.AlignmentEnd < junctionData.Position || readLeftPos - scLength > junctionData.Position)
                return false;

            int juncReadScLength = juncRead.cigar().getFirstCigarElement().getLength();
            int posOffset = juncRead.AlignmentStart - readLeftPos;
            int softClipDiff = juncReadScLength - scLength;
            int junctionReadOffset = softClipDiff - posOffset;
            int juncReadLength = juncRead.readBases().length();

            // check matches from the SC bases up until the end of the first match element or junction/read diff
            int endIndex = scLength + min(max(junctionData.Position - read.AlignmentStart, 0), firstMatchLength);

            int highQualMismatches = 0;
            int baseMatches = 0;

            for(int i = 0; i < endIndex; ++i)
            {
                char readBase = read.readBases().charAt(i);

                int juncIndex = i + junctionReadOffset;
                if(juncIndex < 0 || juncIndex >= juncReadLength)
                    return false;

                char juncReadBase = juncRead.readBases().charAt(juncIndex);

                if(readBase == juncReadBase)
                {
                    ++baseMatches;
                    continue;
                }

                if(belowMinQual(read.baseQualities()[i]) || belowMinQual(juncRead.baseQualities()[juncIndex]))
                    continue;

                ++highQualMismatches;

                if(highQualMismatches > MAX_HIGH_QUAL_BASE_MISMATCHES)
                    return false;
            }

            double baseMatchPerc = baseMatches / (double)endIndex;
            return baseMatchPerc > MIN_EXACT_BASE_PERC;
        }
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

    public static int markSupplementaryDuplicates(final Map<String,ReadGroup> readGroupMap, final ReadIdTrimmer readIdTrimmer)
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
                if(!read.isMateMapped() || !read.hasSuppAlignment())
                    continue;

                if(read.readType() != ReadType.JUNCTION && read.readType() != ReadType.CANDIDATE_SUPPORT)
                    continue;

                int unclippedPosition = unclippedPosition(read);

                List<PrepRead> matchingGroups = initialPositionMap.get(unclippedPosition);

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

                    if(!firstFrag.matches(nextFrag) || firstFrag.IsPrimary == nextFrag.IsPrimary)
                        continue;

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
