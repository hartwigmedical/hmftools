package com.hartwig.hmftools.esvee.prep;

import static java.lang.Math.abs;
import static java.lang.Math.min;

import static com.hartwig.hmftools.esvee.common.CommonUtils.isDiscordantFragment;
import static com.hartwig.hmftools.esvee.prep.PrepConstants.DISCORDANT_GROUP_MIN_FRAGMENTS;
import static com.hartwig.hmftools.esvee.prep.PrepConstants.DISCORDANT_GROUP_MIN_FRAGMENTS_SHORT;
import static com.hartwig.hmftools.esvee.prep.types.DiscordantGroup.firstPrimaryRead;

import java.util.List;
import java.util.Set;

import com.google.common.collect.Lists;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.genome.region.Orientation;
import com.hartwig.hmftools.common.region.ChrBaseRegion;
import com.hartwig.hmftools.esvee.prep.types.DiscordantGroup;
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
                addJunctions(discordantGroup, discordantJunctions);
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
        List<ChrBaseRegion> remoteRegions = discordantGroup.validRemoteRegions(DISCORDANT_GROUP_MIN_FRAGMENTS);

        if(remoteRegions.isEmpty())
            return false;

        for(ChrBaseRegion remoteRegion : remoteRegions)
        {
            if(!discordantGroup.Region.Chromosome.equals(remoteRegion.Chromosome))
                return true;

            int minDistance = min(abs(discordantGroup.innerPosition() - remoteRegion.start()),
                    abs(discordantGroup.innerPosition() - remoteRegion.end()));

            if(minDistance > shortFragmentLength * 2)
                return true;
        }

        return !discordantGroup.validRemoteRegions(DISCORDANT_GROUP_MIN_FRAGMENTS_SHORT).isEmpty();
    }

    private static void addJunctions(final DiscordantGroup discordantGroup, final List<JunctionData> discordantJunctions)
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

        for(ChrBaseRegion remoteRegion : discordantGroup.validRemoteRegions(DISCORDANT_GROUP_MIN_FRAGMENTS))
        {
            junctionData.addRemoteJunction(new RemoteJunction(remoteRegion.Chromosome, remoteRegion.start(), Orientation.FORWARD));
        }
    }

    public static boolean isDiscordantGroup(final ReadGroup readGroup, final int minFragmentLength, final int maxFragmentLength)
    {
        if(readGroup.reads().stream().allMatch(x -> x.record().getSupplementaryAlignmentFlag()))
            return false;

        // only the first read is used and so only that is checked
        return isDiscordantFragment(readGroup.reads().get(0).record(), maxFragmentLength, null);
    }
}
