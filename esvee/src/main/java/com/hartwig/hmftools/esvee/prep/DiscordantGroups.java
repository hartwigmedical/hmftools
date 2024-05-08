package com.hartwig.hmftools.esvee.prep;

import static java.lang.Math.abs;
import static java.lang.Math.min;

import static com.hartwig.hmftools.common.bam.SamRecordUtils.mateNegativeStrand;
import static com.hartwig.hmftools.common.region.BaseRegion.positionWithin;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_END;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_START;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.switchIndex;
import static com.hartwig.hmftools.esvee.prep.PrepConstants.DISCORDANT_GROUP_MAX_DISTANCE;
import static com.hartwig.hmftools.esvee.prep.PrepConstants.DISCORDANT_GROUP_MIN_FRAGMENTS;
import static com.hartwig.hmftools.esvee.prep.PrepConstants.DISCORDANT_GROUP_MIN_FRAGMENTS_SHORT;
import static com.hartwig.hmftools.esvee.prep.types.PrepRead.UNMAPPED_CHR;

import java.util.List;
import java.util.Set;

import com.google.common.collect.Lists;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.region.ChrBaseRegion;
import com.hartwig.hmftools.esvee.prep.types.GroupBoundary;
import com.hartwig.hmftools.esvee.prep.types.JunctionData;
import com.hartwig.hmftools.esvee.prep.types.PrepRead;
import com.hartwig.hmftools.esvee.prep.types.ReadGroup;
import com.hartwig.hmftools.esvee.prep.types.ReadType;
import com.hartwig.hmftools.esvee.prep.types.RemoteJunction;

public final class DiscordantGroups
{
    public static List<JunctionData> formDiscordantJunctions(
            final ChrBaseRegion region, final List<ReadGroup> readGroups, int shortFragmentLength)
    {
        List<JunctionData> discordantJunctions = Lists.newArrayList();
        Set<String> assignedGroups = Sets.newHashSet();

        for(int i = 0; i < readGroups.size() - DISCORDANT_GROUP_MIN_FRAGMENTS;)
        {
            ReadGroup group1 = readGroups.get(i);

            if(assignedGroups.contains(group1.id()))
            {
                ++i;
                continue;
            }

            PrepRead read1 = group1.reads().get(0);

            GroupBoundary[] group1Boundaries = groupBoundaries(group1);
            PrepRead[] boundaryReads = null;
            GroupBoundary[] innerBoundaries = null;
            List<ReadGroup> closeGroups = null;

            int lastSkippedIndex = readGroups.size(); // used to set where the start the next search

            for(int j = i + 1; j < readGroups.size(); ++j)
            {
                ReadGroup group2 = readGroups.get(j);

                if(assignedGroups.contains(group2.id()))
                    continue;

                PrepRead read2 = group2.reads().get(0);

                if(read2.orientation() != read1.orientation() || read2.mateOrientation() != read1.mateOrientation())
                {
                    lastSkippedIndex = min(j, lastSkippedIndex);
                    continue;
                }

                int group2Boundary = read2.orientation().isForward() ? read2.end() : read2.start();

                if(abs(group2Boundary - group1Boundaries[SE_START].Position) > DISCORDANT_GROUP_MAX_DISTANCE)
                {
                    lastSkippedIndex = min(j, lastSkippedIndex);
                    break;
                }

                GroupBoundary[] group2Boundaries = groupBoundaries(group2);

                if(!regionsWithinRange(group1Boundaries, group2Boundaries))
                {
                    lastSkippedIndex = min(j, lastSkippedIndex);
                    continue;
                }

                if(closeGroups == null)
                {
                    closeGroups = Lists.newArrayList(group1);
                    boundaryReads = new PrepRead[] {read1, read1};
                    innerBoundaries = new GroupBoundary[] { group1Boundaries[SE_START], group1Boundaries[SE_END] };
                }

                closeGroups.add(group2);

                // widen with new group and record the reads at the innermost boundary
                if(isCloserToJunction(innerBoundaries, group2Boundaries, SE_START))
                {
                    boundaryReads[SE_START] = read2;
                    innerBoundaries[SE_START] = group2Boundaries[SE_START];
                }

                if(isCloserToJunction(innerBoundaries, group2Boundaries, SE_END))
                {
                    boundaryReads[SE_END] = read2;
                    innerBoundaries[SE_END] = group2Boundaries[SE_END];
                }
            }

            if(closeGroups != null && hasSufficientUnassignedFragments(closeGroups, innerBoundaries, shortFragmentLength))
            {
                addJunctions(closeGroups, innerBoundaries, boundaryReads, region, discordantJunctions);
                closeGroups.forEach(x -> assignedGroups.add(x.id()));

                // jump the first skipped index to avoid reassessing groups already added
                i = lastSkippedIndex == readGroups.size() ? i + 1 : lastSkippedIndex;
            }
            else
            {
                ++i;
            }
        }

        return discordantJunctions;
    }

    private static boolean hasSufficientUnassignedFragments(
            final List<ReadGroup> readGroups, final GroupBoundary[] innerBoundaries, int shortFragmentLength)
    {
        boolean isShortLocalDel = innerBoundaries[SE_START].Chromosome.equals(innerBoundaries[SE_END].Chromosome)
                && innerBoundaries[SE_START].Orient != innerBoundaries[SE_END].Orient
                && abs(innerBoundaries[SE_END].Position - innerBoundaries[SE_START].Position) < shortFragmentLength * 2;

        int minFragments = isShortLocalDel ? DISCORDANT_GROUP_MIN_FRAGMENTS_SHORT : DISCORDANT_GROUP_MIN_FRAGMENTS;
        return readGroups.size() >= minFragments && readGroups.stream().filter(x -> !x.hasJunctionPositions()).count() >= minFragments;
    }

    private static void addJunctions(
            final List<ReadGroup> readGroups, final GroupBoundary[] innerBoundaries, final PrepRead[] boundaryReads,
            final ChrBaseRegion region, final List<JunctionData> discordantJunctions)
    {
        for(int se = SE_START; se <= SE_END; ++se)
        {
            // only create junctions from groups within this region
            if(!region.containsPosition(innerBoundaries[se].Chromosome, innerBoundaries[se].Position))
                continue;

            JunctionData junctionData = new JunctionData(innerBoundaries[se].Position, innerBoundaries[se].Orient, boundaryReads[se]);
            discordantJunctions.add(junctionData);

            junctionData.markDiscordantGroup();

            for(ReadGroup readGroup : readGroups)
            {
                junctionData.SupportingGroups.add(readGroup);
                readGroup.addJunctionPosition(junctionData);

                readGroup.reads().forEach(x -> x.setReadType(ReadType.SUPPORT, true));
                // readGroup.reads().forEach(x -> junctionData.addReadType(x, ReadType.SUPPORT)); // no need
            }

            int seOther = switchIndex(se);

            junctionData.addRemoteJunction(new RemoteJunction(
                    innerBoundaries[seOther].Chromosome, innerBoundaries[seOther].Position, innerBoundaries[seOther].Orient));
        }
    }

    public static boolean isDiscordantGroup(final ReadGroup readGroup, final int minFragmentLength, final int maxFragmentLength)
    {
        /* added to support unpaired reads but prevents a single read spanning partitions from forming a discordant group
        if(readGroup.reads().stream().filter(x -> !x.isSupplementaryAlignment()).count() < 2)
            return false;
        */

        // only the first read is used and so only that is checked
        return isDiscordantRead(readGroup.reads().get(0), minFragmentLength, maxFragmentLength);
    }

    private static boolean isDiscordantRead(final PrepRead read, final int minFragmentLength, final int maxFragmentLength)
    {
        if(read.Chromosome.equals(UNMAPPED_CHR) || read.MateChromosome.equals(UNMAPPED_CHR))
            return false;

        if(read.fragmentInsertSize() < minFragmentLength || read.fragmentInsertSize() > maxFragmentLength)
            return true;

        if(!read.Chromosome.equals(read.MateChromosome))
            return true;

        if(read.record().getReadNegativeStrandFlag() == mateNegativeStrand(read.record()))
            return true;

        return false;
    }

    private static boolean isCloserToJunction(final GroupBoundary[] current, final GroupBoundary[] test, int seIndex)
    {
        if(current[seIndex].Orient.isForward())
        {
            return test[seIndex].Position > current[seIndex].Position;
        }
        else
        {
            return test[seIndex].Position < current[seIndex].Position;
        }
    }

    private static boolean regionsWithinRange(final GroupBoundary[] boundaries1, final GroupBoundary[] boundaries2)
    {
        for(int se = SE_START; se <= SE_END; ++se)
        {
            if(!boundaries1[se].Chromosome.equals(boundaries2[se].Chromosome))
                return false;

            if(!positionWithin(
                    boundaries2[se].Position,
                    boundaries1[se].Position - DISCORDANT_GROUP_MAX_DISTANCE,
                    boundaries1[se].Position + DISCORDANT_GROUP_MAX_DISTANCE))
            {
                return false;
            }
        }

        return true;
    }

    private static GroupBoundary[] groupBoundaries(final ReadGroup readGroup)
    {
        PrepRead read = readGroup.reads().get(0);

        GroupBoundary boundary1 = new GroupBoundary(
                read.Chromosome, read.orientation().isForward() ? read.end() : read.start(),
                read.orientation());

        GroupBoundary boundary2;

        if(readGroup.size() == 2)
        {
            PrepRead read2 = readGroup.reads().get(1);
            boundary2 = new GroupBoundary(
                    read2.Chromosome, read2.orientation().isForward() ? read2.end() : read2.start(),
                    read2.orientation());
        }
        else
        {
            boundary2 = new GroupBoundary(
                    read.MateChromosome,
                    read.mateOrientation().isForward() ? read.MatePosStart + read.record().getReadLength() : read.MatePosStart,
                    read.mateOrientation());
        }

        return new GroupBoundary[] { boundary1, boundary2 };
    }
}
