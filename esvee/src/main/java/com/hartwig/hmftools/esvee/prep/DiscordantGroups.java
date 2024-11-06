package com.hartwig.hmftools.esvee.prep;

import static java.lang.Math.abs;
import static java.lang.Math.min;

import static com.hartwig.hmftools.esvee.common.CommonUtils.isDiscordantFragment;
import static com.hartwig.hmftools.esvee.common.SvConstants.MIN_MAP_QUALITY;
import static com.hartwig.hmftools.esvee.prep.KnownHotspot.readGroupMatchesHotspot;
import static com.hartwig.hmftools.esvee.prep.PrepConstants.DISCORDANT_GROUP_MAX_LOCAL_LENGTH;
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
import java.util.stream.Collectors;

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

public class DiscordantGroups
{
    private final ChrBaseRegion mRegion;
    private final int mMaxConcordantFragmentLength;
    private final boolean mTrackRemotes;
    private final List<KnownHotspot> mKnownHotspots;

    public DiscordantGroups(final ChrBaseRegion region, int shortFragmentLength, final List<KnownHotspot> knownHotspots, boolean trackRemotes)
    {
        mRegion = region;
        mMaxConcordantFragmentLength = shortFragmentLength;
        mKnownHotspots = knownHotspots;
        mTrackRemotes = trackRemotes;
    }

    public List<JunctionData> formDiscordantJunctions(final List<ReadGroup> candidateReadGroups)
    {
        List<JunctionData> discordantJunctions = Lists.newArrayList();

        // each starting orientation is tested in turn since these need to be consistent to form a candidate discordant-only junction
        for(int o = 0; o <= 1; ++o)
        {
            Orientation orientation = o == 0 ? Orientation.FORWARD : Orientation.REVERSE;

            List<ReadGroup> readGroups = candidateReadGroups.stream()
                    .filter(x -> firstPrimaryRead(x).orientation() == orientation).collect(Collectors.toList());

            Set<String> assignedGroups = Sets.newHashSet();

            // first sort the groups by their first primary read coordinates, which will be used to form discordant group boundaries
            Collections.sort(readGroups, new ReadGroupSorter());

            for(int i = 0; i < readGroups.size() - DISCORDANT_GROUP_MIN_FRAGMENTS; )
            {
                ReadGroup firstGroup = readGroups.get(i);

                if(assignedGroups.contains(firstGroup.id()))
                {
                    ++i;
                    continue;
                }

                DiscordantGroup discordantGroup = DiscordantGroup.fromReadGroup(firstGroup);

                int j = i + 1;
                for(; j < readGroups.size(); ++j)
                {
                    ReadGroup nextGroup = readGroups.get(j);

                    if(assignedGroups.contains(nextGroup.id()))
                        continue;

                    PrepRead nextRead = firstPrimaryRead(nextGroup);

                    if(nextRead.start() > discordantGroup.Region.end())
                        break;

                    discordantGroup.tryAddReadGroup(nextGroup, nextRead);
                }

                if(isValidDiscordantGroup(discordantGroup))
                {
                    addJunctions(discordantGroup, discordantJunctions);
                    discordantGroup.readGroups().forEach(x -> assignedGroups.add(x.id()));
                }

                i = j; // move to the first group not added
            }
        }

        return discordantJunctions;
    }

    private boolean isValidDiscordantGroup(final DiscordantGroup discordantGroup)
    {
        if(discordantGroup.readGroups().size() < DISCORDANT_GROUP_MIN_FRAGMENTS)
            return false;

        // only consider junctions from groups within this region
        if(!mRegion.overlaps(discordantGroup.Region))
            return false;

        // only consider pairings which are either local within the defined distance limit or in a known fusion region

        // check that at least one remote region also has sufficient reads
        List<DiscordantRemoteRegion> remoteRegions = discordantGroup.remoteRegions();

        Collections.sort(remoteRegions, Comparator.comparingInt(x -> -x.readCount()));

        if(remoteRegions.get(0).readCount() < DISCORDANT_GROUP_MIN_FRAGMENTS)
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
            if(hasRemoteRequiredFragments(discordantGroup, remoteRegion))
                return true;
        }

        return false;
    }

    private void addJunctions(final DiscordantGroup discordantGroup, final List<JunctionData> discordantJunctions)
    {
        // define the junction point from the inner most read from the remote region with the most support
        DiscordantRemoteRegion mainRemoteRegion = discordantGroup.remoteRegions().get(0);

        PrepRead innerRead = null;
        int innerPosition = 0;

        for(ReadGroup readGroup : mainRemoteRegion.ReadGroups)
        {
            PrepRead read = firstPrimaryRead(readGroup);

            if(discordantGroup.Orient.isForward())
            {
                if(innerRead == null || read.end() > innerRead.end())
                {
                    innerRead = read;
                    innerPosition = read.end();
                }
            }
            else
            {
                if(innerRead == null || read.start() < innerRead.start())
                {
                    innerRead = read;
                    innerPosition = read.start();
                }
            }
        }

        JunctionData junctionData = new JunctionData(innerPosition, discordantGroup.Orient, innerRead);
        discordantGroup.setInnerRead(innerRead);

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
            if(mTrackRemotes)
            {
                RemoteJunction remoteJunction = new RemoteJunction(remoteRegion.Chromosome, remoteRegion.start(), Orientation.FORWARD);
                remoteJunction.Fragments = remoteRegion.readCount();
                junctionData.addRemoteJunction(remoteJunction);
            }

            if(!hasRemoteRequiredFragments(discordantGroup, remoteRegion))
                continue;

            int remoteJunctionPosition = remoteRegion.Orient.isForward() ? remoteRegion.end() : remoteRegion.start();

            if(!mRegion.containsPosition(remoteRegion.Chromosome, remoteJunctionPosition))
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

            if(mTrackRemotes)
            {
                RemoteJunction remoteJunction = new RemoteJunction(
                        discordantGroup.Region.Chromosome, discordantGroup.innerPosition(), discordantGroup.Orient);
                remoteJunction.Fragments = discordantGroup.readGroups().size();
                remoteJunctionData.addRemoteJunction(remoteJunction);
            }
        }
    }

    private boolean hasRemoteRequiredFragments(final DiscordantGroup discordantGroup, final DiscordantRemoteRegion remoteRegion)
    {
        boolean isShortLocal = false;

        if(discordantGroup.Region.Chromosome.equals(remoteRegion.Chromosome))
        {
            int minDistance = remoteRegion.start() > discordantGroup.Region.end() ?
                    remoteRegion.start() - discordantGroup.Region.end() : discordantGroup.Region.start() - remoteRegion.end();

            isShortLocal = minDistance <= mMaxConcordantFragmentLength * 2;
        }

        if(isShortLocal)
            return remoteRegion.readCount() >= DISCORDANT_GROUP_MIN_FRAGMENTS_SHORT;
        else
            return remoteRegion.readCount() >= DISCORDANT_GROUP_MIN_FRAGMENTS;
    }

    public boolean isDiscordantGroup(final ReadGroup readGroup)
    {
        boolean hasNonSupp = false;
        boolean hasPassingMapQual = false;

        for(PrepRead read : readGroup.reads())
        {
            hasNonSupp |= !read.isSupplementaryAlignment();
            hasPassingMapQual |= read.record().getMappingQuality() >= MIN_MAP_QUALITY;
        }

        if(!hasNonSupp || !hasPassingMapQual)
            return false;

        // only the first read is used and so only that is checked
        PrepRead firstRead = readGroup.reads().stream().filter(x -> !x.isSupplementaryAlignment()).findFirst().orElse(null);

        if(firstRead == null)
            return false;

        if(!isDiscordantFragment(firstRead.record(), mMaxConcordantFragmentLength, null))
            return false;

        // must then either be in a known hotspot pair or a local DEL or DUP within the required distance
        if(readGroupMatchesHotspot(mKnownHotspots, readGroup))
            return true;

        // must be local and a DEL or DUP
        if(!firstRead.Chromosome.equals(firstRead.MateChromosome) || firstRead.orientation() == firstRead.mateOrientation())
            return false;

        int length = abs(firstRead.start() - firstRead.record().getMateAlignmentStart());
        return length <= DISCORDANT_GROUP_MAX_LOCAL_LENGTH;
    }

    private static class ReadGroupSorter implements Comparator<ReadGroup>
    {
        public int compare(final ReadGroup first, final ReadGroup second)
        {
            PrepRead firstRead = firstPrimaryRead(first);
            PrepRead secondRead = firstPrimaryRead(second);
            return Integer.compare(firstRead.start(), secondRead.start());
        }
    }

}
