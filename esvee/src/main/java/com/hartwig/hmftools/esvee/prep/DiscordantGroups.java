package com.hartwig.hmftools.esvee.prep;

import static java.lang.Math.abs;
import static java.lang.Math.max;
import static java.lang.Math.min;

import static com.hartwig.hmftools.esvee.common.CommonUtils.isDiscordantFragment;
import static com.hartwig.hmftools.esvee.prep.PrepConstants.DISCORDANT_GROUP_MIN_ALIGN_SCORE;
import static com.hartwig.hmftools.esvee.prep.PrepConstants.DISCORDANT_GROUP_MIN_FRAGMENTS;
import static com.hartwig.hmftools.esvee.prep.PrepConstants.DISCORDANT_GROUP_MIN_FRAGMENTS_SHORT;
import static com.hartwig.hmftools.esvee.prep.PrepConstants.DISCORDANT_GROUP_MIN_MAP_QUAL;
import static com.hartwig.hmftools.esvee.prep.ReadFilters.aboveRepeatTrimmedAlignmentThreshold;
import static com.hartwig.hmftools.esvee.prep.types.DiscordantGroup.firstPrimaryRead;

import java.util.Collections;
import java.util.Comparator;
import java.util.List;
import java.util.Set;

import com.google.common.collect.Lists;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.genome.region.Orientation;
import com.hartwig.hmftools.common.region.ChrBaseRegion;
import com.hartwig.hmftools.esvee.prep.types.DiscordantGroup;
import com.hartwig.hmftools.esvee.prep.types.DiscordantRemoteRegion;
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

        // first sort the groups by their first primary read coordinates, which will be used to form discordant group boundaries
        Collections.sort(readGroups, new ReadGroupSorter());

        for(int i = 0; i < readGroups.size() - DISCORDANT_GROUP_MIN_FRAGMENTS;)
        {
            ReadGroup firstGroup = readGroups.get(i);

            if(assignedGroups.contains(firstGroup.id()) || firstGroup.hasJunctionPositions())
            {
                ++i;
                continue;
            }

            DiscordantGroup discordantGroup = DiscordantGroup.fromReadGroup(firstGroup);

            int lastSkippedIndex = readGroups.size(); // used to set where the start the next search

            for(int j = i + 1; j < readGroups.size(); ++j)
            {
                ReadGroup nextGroup = readGroups.get(j);

                if(assignedGroups.contains(nextGroup.id()))
                    continue;

                PrepRead nextRead = firstPrimaryRead(nextGroup);

                if(!discordantGroup.tryAddReadGroup(nextGroup, nextRead))
                {
                    lastSkippedIndex = min(j, lastSkippedIndex);

                    if(nextRead.start() > discordantGroup.Region.end())
                        break;
                }
            }

            if(isValidDiscordantGroup(region, discordantGroup, shortFragmentLength))
            {
                addJunctions(region, discordantGroup, discordantJunctions, shortFragmentLength);
                discordantGroup.readGroups().forEach(x -> assignedGroups.add(x.id()));

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

    private static boolean isValidDiscordantGroup(final ChrBaseRegion region, final DiscordantGroup discordantGroup, int shortFragmentLength)
    {
        // only consider junctions from groups within this region
        if(!region.containsPosition(discordantGroup.Region.Chromosome, discordantGroup.innerPosition()))
            return false;

        if(discordantGroup.readGroups().size() < DISCORDANT_GROUP_MIN_FRAGMENTS)
            return false;

        // check that at least one remote region also has sufficient reads
        List<DiscordantRemoteRegion> remoteRegions = discordantGroup.remoteRegions();

        if(remoteRegions.stream().noneMatch(x -> x.readCount() >= DISCORDANT_GROUP_MIN_FRAGMENTS))
            return false;

        boolean aboveMinMapQual = false;
        boolean aboveMinAlignScore = false;

        for(ReadGroup readGroup : discordantGroup.readGroups())
        {
            PrepRead read = firstPrimaryRead(readGroup);

            aboveMinMapQual |= read.mapQuality() >= DISCORDANT_GROUP_MIN_MAP_QUAL;

            aboveMinAlignScore |= aboveRepeatTrimmedAlignmentThreshold(read, DISCORDANT_GROUP_MIN_ALIGN_SCORE, false);

            if(aboveMinAlignScore && aboveMinMapQual)
                break;
        }

        if(!aboveMinAlignScore || !aboveMinMapQual)
            return false;

        for(DiscordantRemoteRegion remoteRegion : remoteRegions)
        {
            if(hasRemoteRequiredFragments(discordantGroup, remoteRegion, shortFragmentLength))
                return true;
        }

        return false;
    }

    private static void addJunctions(
            final ChrBaseRegion region, final DiscordantGroup discordantGroup, final List<JunctionData> discordantJunctions,
            int shortFragmentLength)
    {
        JunctionData junctionData = new JunctionData(
                discordantGroup.innerPosition(), discordantGroup.Orient, discordantGroup.innerRead());
        discordantJunctions.add(junctionData);

        junctionData.markDiscordantGroup();

        for(ReadGroup readGroup : discordantGroup.readGroups())
        {
            junctionData.SupportingGroups.add(readGroup);
            readGroup.addJunctionPosition(junctionData);

            readGroup.reads().forEach(x -> x.setReadType(ReadType.SUPPORT, true));
        }

        List<DiscordantRemoteRegion> remoteRegions = discordantGroup.remoteRegions();

        // create junctions from remote regions which satisfy the required fragment count
        for(DiscordantRemoteRegion remoteRegion : remoteRegions)
        {
            junctionData.addRemoteJunction(new RemoteJunction(remoteRegion.Chromosome, remoteRegion.start(), Orientation.FORWARD));

            if(!hasRemoteRequiredFragments(discordantGroup, remoteRegion, shortFragmentLength))
                continue;

            int remoteJunctionPosition = remoteRegion.Orient.isForward() ? remoteRegion.end() : remoteRegion.start();

            if(!region.containsPosition(remoteRegion.Chromosome, remoteJunctionPosition))
                continue;

            JunctionData remoteJunctionData = new JunctionData(
                    remoteJunctionPosition, remoteRegion.Orient, remoteRegion.ReadGroups.get(0).reads().get(0));

            discordantJunctions.add(remoteJunctionData);

            remoteJunctionData.markDiscordantGroup();

            for(ReadGroup readGroup : remoteRegion.ReadGroups)
            {
                remoteJunctionData.SupportingGroups.add(readGroup);
                readGroup.addJunctionPosition(remoteJunctionData);

                readGroup.reads().forEach(x -> x.setReadType(ReadType.SUPPORT, true));
            }
        }
    }

    private static boolean hasRemoteRequiredFragments(
            final DiscordantGroup discordantGroup, final DiscordantRemoteRegion remoteRegion, int shortFragmentLength)
    {
        boolean isShortLocal = discordantGroup.Region.Chromosome.equals(remoteRegion.Chromosome)
                && min(abs(discordantGroup.innerPosition() - remoteRegion.start()),
                abs(discordantGroup.innerPosition() - remoteRegion.end())) <= shortFragmentLength * 2;

        if(isShortLocal)
            return remoteRegion.readCount() >= DISCORDANT_GROUP_MIN_FRAGMENTS_SHORT;
        else
            return remoteRegion.readCount() >= DISCORDANT_GROUP_MIN_FRAGMENTS;
    }

    public static boolean isDiscordantGroup(final ReadGroup readGroup, final int maxFragmentLength)
    {
        if(readGroup.reads().stream().allMatch(x -> x.record().getSupplementaryAlignmentFlag()))
            return false;

        // only the first read is used and so only that is checked
        return isDiscordantFragment(readGroup.reads().get(0).record(), maxFragmentLength, null);
    }

    public static class ReadGroupSorter implements Comparator<ReadGroup>
    {
        public int compare(final ReadGroup first, final ReadGroup second)
        {
            PrepRead firstRead = firstPrimaryRead(first);
            PrepRead secondRead = firstPrimaryRead(second);
            return Integer.compare(firstRead.start(), secondRead.start());
        }
    }

}
