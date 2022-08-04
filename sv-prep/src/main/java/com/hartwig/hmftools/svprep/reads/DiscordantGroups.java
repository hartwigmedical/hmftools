package com.hartwig.hmftools.svprep.reads;

import static java.lang.Math.abs;

import static com.hartwig.hmftools.common.utils.sv.BaseRegion.positionWithin;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_END;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_START;
import static com.hartwig.hmftools.common.utils.sv.SvCommonUtils.POS_ORIENT;
import static com.hartwig.hmftools.svprep.reads.ReadRecord.UNMAPPED_CHR;

import java.util.List;
import java.util.Set;

import com.google.common.collect.Lists;
import com.google.common.collect.Sets;

public final class DiscordantGroups
{
    public static final int MIN_FRAGMENT_COUNT = 5;
    public static final int MAX_START_DISTANCE = 500;
    public static final int MAX_END_DISTANCE = 1000;

    public static List<JunctionData> formDiscordantJunctions(final List<ReadGroup> readGroups)
    {
        List<JunctionData> discordantJunctions = Lists.newArrayList();
        Set<String> assignedGroups = Sets.newHashSet();

        for(int i = 0; i < readGroups.size() - MIN_FRAGMENT_COUNT; ++i)
        {
            ReadGroup group1 = readGroups.get(i);

            GroupBoundary[] group1Boundaries = groupBoundaries(group1);
            ReadRecord read1 = group1.reads().get(0);
            ReadRecord[] boundaryReads = null;
            GroupBoundary[] innerBoundaries = null;
            List<ReadGroup> closeGroups = null;

            for(int j = i + 1; j < readGroups.size(); ++j)
            {
                ReadGroup group2 = readGroups.get(j);
                ReadRecord read2 = group2.reads().get(0);

                if(read2.orientation() != read1.orientation() || read2.mateOrientation() != read1.mateOrientation())
                    continue;

                int group2Boundary = read2.orientation() == POS_ORIENT ? read2.end() : read2.start();

                if(abs(group2Boundary - group1Boundaries[SE_START].Position) > MAX_START_DISTANCE)
                    break;

                if(assignedGroups.contains(group2.id()))
                    continue;

                GroupBoundary[] group2Boundaries = groupBoundaries(group2);

                if(!regionsWithinRange(group1Boundaries, group2Boundaries))
                    continue;

                if(closeGroups == null)
                {
                    closeGroups = Lists.newArrayList(group1);
                    boundaryReads = new ReadRecord[] {read1, read1};
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

            if(closeGroups != null && closeGroups.size() >= MIN_FRAGMENT_COUNT)
            {
                closeGroups.forEach(x -> assignedGroups.add(x.id()));
                addJunctions(closeGroups, innerBoundaries, boundaryReads, discordantJunctions);
            }
        }

        return discordantJunctions;
    }

    private static void addJunctions(
            final List<ReadGroup> readGroups, final GroupBoundary[] innerBoundaries, final ReadRecord[] boundaryReads,
            final List<JunctionData> discordantJunctions)
    {
        for(int se = SE_START; se <= SE_END; ++se)
        {
            JunctionData junctionData = new JunctionData(innerBoundaries[se].Position, innerBoundaries[se].Orientation, boundaryReads[se]);
            discordantJunctions.add(junctionData);

            junctionData.markDiscordantGroup();

            for(ReadGroup readGroup : readGroups)
            {
                junctionData.SupportingGroups.add(readGroup);
                readGroup.addJunctionPosition(junctionData.Position);

                readGroup.reads().forEach(x -> x.setReadType(ReadType.SUPPORT, true));
                // readGroup.reads().forEach(x -> junctionData.addReadType(x, ReadType.SUPPORT)); // no need
            }
        }
    }

    public static boolean isDiscordantGroup(final ReadGroup readGroup, final int maxFragmentLength)
    {
        // only the first read is used and so only that is checked
        return isDiscordantRead(readGroup.reads().get(0), maxFragmentLength);
    }

    private static boolean isDiscordantRead(final ReadRecord read, final int maxFragmentLength)
    {
        if(read.Chromosome.equals(UNMAPPED_CHR) || read.MateChromosome.equals(UNMAPPED_CHR))
            return false;

        if(read.fragmentInsertSize() > maxFragmentLength)
            return true;

        if(!read.Chromosome.equals(read.MateChromosome))
            return true;

        if(read.record().getReadNegativeStrandFlag() == read.record().getMateNegativeStrandFlag())
            return true;

        return false;
    }

    private static boolean isCloserToJunction(final GroupBoundary[] current, final GroupBoundary[] test, int seIndex)
    {
        if(current[seIndex].Orientation == POS_ORIENT)
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
                    boundaries1[se].Position - MAX_START_DISTANCE,
                    boundaries1[se].Position + MAX_START_DISTANCE))
            {
                return false;
            }
        }

        return true;
    }

    private static GroupBoundary[] groupBoundaries(final ReadGroup readGroup)
    {
        ReadRecord read = readGroup.reads().get(0);

        GroupBoundary boundary1 = new GroupBoundary(
                read.Chromosome, read.orientation() == POS_ORIENT ? read.end() : read.start(),
                read.orientation());

        GroupBoundary boundary2;

        if(readGroup.size() == 2)
        {
            ReadRecord read2 = readGroup.reads().get(1);
            boundary2 = new GroupBoundary(
                    read2.Chromosome, read2.orientation() == POS_ORIENT ? read2.end() : read2.start(),
                    read2.orientation());
        }
        else
        {
            boundary2 = new GroupBoundary(
                    read.MateChromosome,
                    read.mateOrientation() == POS_ORIENT ? read.MatePosStart + read.record().getReadLength() : read.MatePosStart,
                    read.mateOrientation());
        }

        return new GroupBoundary[] { boundary1, boundary2 };
    }
}
