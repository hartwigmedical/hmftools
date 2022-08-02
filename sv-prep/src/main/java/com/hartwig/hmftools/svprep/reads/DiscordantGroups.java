package com.hartwig.hmftools.svprep.reads;

import static java.lang.Math.max;
import static java.lang.Math.min;

import static com.hartwig.hmftools.common.utils.sv.BaseRegion.positionWithin;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_END;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_START;
import static com.hartwig.hmftools.common.utils.sv.SvCommonUtils.NEG_ORIENT;
import static com.hartwig.hmftools.common.utils.sv.SvCommonUtils.POS_ORIENT;
import static com.hartwig.hmftools.svprep.SvCommon.SV_LOGGER;
import static com.hartwig.hmftools.svprep.reads.ReadRecord.UNMAPPED_CHR;

import java.util.List;
import java.util.Set;

import com.beust.jcommander.internal.Sets;
import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.utils.sv.ChrBaseRegion;

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

            ChrBaseRegion[] regions1 = groupRegions(group1);
            List<ReadGroup> closeGroups = null;

            for(int j = i + 1; j < readGroups.size(); ++j)
            {
                ReadGroup group2 = readGroups.get(j);

                if(assignedGroups.contains(group2.id()))
                    continue;

                ChrBaseRegion[] regions2 = groupRegions(group2);

                if(regionsWithinRange(regions1, regions2))
                {
                    if(closeGroups == null)
                        closeGroups = Lists.newArrayList(group1);

                    closeGroups.add(group2);

                    // widen with new group
                    regions1[SE_START].setStart(min(regions1[SE_START].start(), regions2[SE_START].start()));
                    regions1[SE_START].setEnd(max(regions1[SE_START].end(), regions2[SE_START].end()));
                    regions1[SE_END].setStart(min(regions1[SE_END].start(), regions2[SE_END].start()));
                    regions1[SE_END].setEnd(max(regions1[SE_END].end(), regions2[SE_END].end()));
                }
            }

            if(closeGroups != null && closeGroups.size() >= MIN_FRAGMENT_COUNT)
            {
                closeGroups.forEach(x -> assignedGroups.add(x.id()));
                addJunctions(closeGroups, regions1, discordantJunctions);
            }
        }

        return discordantJunctions;
    }

    private static void addJunctions(
            final List<ReadGroup> readGroups, final ChrBaseRegion[] regions, final List<JunctionData> discordantJunctions)
    {
        // determine orientation and then create junctions for any local regions
        for(int se = SE_START; se <= SE_END; ++se)
        {
            // find a read matching the region boundary and orientation
            ChrBaseRegion region = regions[se];
            int junctionPosition = 0;
            byte junctionOrientation = 0;
            ReadRecord boundaryRead = null;

            for(ReadGroup readGroup : readGroups)
            {
                ReadRecord read = readGroup.reads().stream()
                        .filter(x -> x.start() == region.start() && x.orientation() == NEG_ORIENT).findFirst().orElse(null);

                if(read != null)
                {
                    junctionPosition = read.start();
                    junctionOrientation = read.orientation();
                    boundaryRead = read;
                    break;
                }

                read = readGroup.reads().stream()
                        .filter(x -> x.start() == region.end() && x.orientation() == POS_ORIENT).findFirst().orElse(null);

                if(read != null)
                {
                    junctionPosition = read.end();
                    junctionOrientation = read.orientation();
                    boundaryRead = read;
                    break;
                }
            }

            if(boundaryRead != null)
            {
                JunctionData junctionData = new JunctionData(junctionPosition, junctionOrientation, boundaryRead);
                discordantJunctions.add(junctionData);

                junctionData.markDiscordantGroup();
                readGroups.forEach(x -> junctionData.SupportingGroups.add(x));
                readGroups.forEach(x -> x.addJunctionPosition(junctionData.Position));
            }
        }
    }

    public static boolean isDiscordantGroup(final ReadGroup readGroup, final int maxFragmentLength)
    {
        // only the first read is used and so only that is checked
        return isDiscordantRead(readGroup.reads().get(0), maxFragmentLength);
        // return readGroup.reads().stream().anyMatch(x -> isDiscordantRead(x, maxFragmentLength));
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

    private static boolean regionsWithinRange(final ChrBaseRegion[] regions1, final ChrBaseRegion[] regions2)
    {
        for(int se = SE_START; se <= SE_END; ++se)
        {
            if(!regions1[se].Chromosome.equals(regions2[se].Chromosome))
                return false;

            if(!positionWithin(
                    regions2[se].start(),
                    regions1[se].start() - MAX_START_DISTANCE,
                    regions1[se].start() + MAX_START_DISTANCE))
            {
                return false;
            }
        }

        return true;
    }

    private static ChrBaseRegion[] groupRegions(final ReadGroup readGroup)
    {
        ReadRecord read = readGroup.reads().get(0);
        String chr1 = read.Chromosome;
        String chr2 = read.MateChromosome;
        int pos1 = read.orientation() == POS_ORIENT ? read.end() : read.start();
        int pos2 = readGroup.reads().get(0).MatePosStart;


        boolean firstIsLower = false;

        if(HumanChromosome.chromosomeRank(chr1) < HumanChromosome.chromosomeRank(chr2))
        {
            firstIsLower = true;
        }
        else if(HumanChromosome.chromosomeRank(chr1) > HumanChromosome.chromosomeRank(chr2))
        {
            firstIsLower = false;
        }
        else
        {
            firstIsLower = pos1 < pos2;
        }

        return firstIsLower ?
                new ChrBaseRegion[] { new ChrBaseRegion(chr1, pos1, pos1), new ChrBaseRegion(chr2, pos2, pos2) } :
                new ChrBaseRegion[] { new ChrBaseRegion(chr2, pos2, pos2), new ChrBaseRegion(chr1, pos1, pos1) };
    }
}
