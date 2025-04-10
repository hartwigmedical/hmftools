package com.hartwig.hmftools.esvee.prep;

import static java.lang.Math.abs;
import static java.lang.Math.max;
import static java.lang.Math.min;

import static com.hartwig.hmftools.common.region.BaseRegion.positionWithin;
import static com.hartwig.hmftools.esvee.assembly.AssemblyConfig.SV_LOGGER;
import static com.hartwig.hmftools.esvee.common.CommonUtils.belowMinQual;
import static com.hartwig.hmftools.esvee.prep.PrepConstants.MAX_HIGH_QUAL_BASE_MISMATCHES;
import static com.hartwig.hmftools.esvee.prep.PrepConstants.MIN_EXACT_BASE_PERC;
import static com.hartwig.hmftools.esvee.prep.ReadFilters.isChimericRead;
import static com.hartwig.hmftools.esvee.prep.types.FragmentData.unclippedPosition;
import static com.hartwig.hmftools.esvee.prep.types.ReadGroupStatus.DUPLICATE;

import static htsjdk.samtools.CigarOperator.M;
import static htsjdk.samtools.CigarOperator.S;

import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
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
    public static boolean hasOtherJunctionSupport(
            final PrepRead read, final JunctionData junctionData, final ReadFilterConfig filterConfig)
    {
        int unclippedStart = read.unclippedStart();
        int unclippedEnd = read.unclippedEnd();

        // first check for a read crossing the junction
        if(positionWithin(junctionData.Position, unclippedStart, unclippedEnd))
        {
            // correct side of the junction
            int junctionDistance = 0;

            if(junctionData.isForward())
            {
                junctionDistance = min(abs(unclippedEnd - junctionData.Position), abs(read.end() - junctionData.Position));
            }
            else
            {
                junctionDistance = min(abs(unclippedStart - junctionData.Position), abs(read.start() - junctionData.Position));
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
            if(read.end() > junctionData.Position)
                return false;

            junctionDistance = abs(read.end() - junctionData.Position);
        }
        else
        {
            if(read.start() < junctionData.Position) //  || abs(read.end() - junctionData.Position) > filterConfig.maxSupportingFragmentDistance()
                return false;

            junctionDistance = abs(read.start() - junctionData.Position);
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

            int readRightPos = read.end();

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
            if(read.start() > junctionData.Position || readRightPos + scLength < junctionData.Position)
                return false;

            int readEndPosIndex = readLength - scLength - 1;

            int juncReadLength = juncRead.readBases().length();
            int juncReadScLength = juncRead.cigar().getLastCigarElement().getLength();
            int juncReadEndPosIndex = juncReadLength - juncReadScLength - 1;
            int endPosDiff = juncRead.end() - readRightPos;

            int junctionReadOffset = juncReadEndPosIndex - readEndPosIndex - endPosDiff;

            // test all overlapping bases - either from ref or soft-clip bases
            int startIndex = readLength - scLength - min(max(read.end() - junctionData.Position, 0), firstMatchLength);

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

            int readLeftPos = read.start();

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

            if(read.end() < junctionData.Position || readLeftPos - scLength > junctionData.Position)
                return false;

            int juncReadScLength = juncRead.cigar().getFirstCigarElement().getLength();
            int posOffset = juncRead.start() - readLeftPos;
            int softClipDiff = juncReadScLength - scLength;
            int junctionReadOffset = softClipDiff - posOffset;
            int juncReadLength = juncRead.readBases().length();

            // check matches from the SC bases up until the end of the first match element or junction/read diff
            int endIndex = scLength + min(max(junctionData.Position - read.start(), 0), firstMatchLength);

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
