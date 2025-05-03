package com.hartwig.hmftools.esvee.prep;

import static java.lang.Math.abs;
import static java.lang.Math.max;
import static java.lang.Math.min;

import static com.hartwig.hmftools.common.genome.region.Orientation.FORWARD;
import static com.hartwig.hmftools.common.genome.region.Orientation.REVERSE;
import static com.hartwig.hmftools.common.region.BaseRegion.positionsOverlap;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_END;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_START;
import static com.hartwig.hmftools.esvee.assembly.AssemblyConfig.SV_LOGGER;
import static com.hartwig.hmftools.esvee.common.SvConstants.LOW_BASE_QUAL_THRESHOLD;
import static com.hartwig.hmftools.esvee.prep.JunctionUtils.hasExactJunctionSupport;
import static com.hartwig.hmftools.esvee.prep.JunctionUtils.hasOtherJunctionSupport;
import static com.hartwig.hmftools.esvee.prep.JunctionUtils.markSupplementaryDuplicates;
import static com.hartwig.hmftools.esvee.prep.KnownHotspot.junctionMatchesHotspot;
import static com.hartwig.hmftools.esvee.prep.PrepConstants.MIN_HOTSPOT_JUNCTION_SUPPORT;
import static com.hartwig.hmftools.esvee.prep.PrepConstants.MIN_LINE_SOFT_CLIP_LENGTH;
import static com.hartwig.hmftools.esvee.prep.PrepConstants.UNPAIRED_READ_JUNCTION_DISTANCE;
import static com.hartwig.hmftools.esvee.prep.types.ReadFilterType.INSERT_MAP_OVERLAP;
import static com.hartwig.hmftools.esvee.prep.types.ReadFilterType.SOFT_CLIP_LENGTH;
import static com.hartwig.hmftools.esvee.prep.ReadFilters.aboveRepeatTrimmedAlignmentThreshold;
import static com.hartwig.hmftools.esvee.prep.types.ReadType.NO_SUPPORT;

import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.stream.Collectors;

import com.beust.jcommander.internal.Sets;
import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.genome.region.Orientation;
import com.hartwig.hmftools.common.perf.PerformanceCounter;
import com.hartwig.hmftools.common.region.BaseRegion;
import com.hartwig.hmftools.common.region.ChrBaseRegion;
import com.hartwig.hmftools.esvee.prep.types.DiscordantStats;
import com.hartwig.hmftools.esvee.prep.types.JunctionData;
import com.hartwig.hmftools.esvee.prep.types.ReadFilterConfig;
import com.hartwig.hmftools.esvee.prep.types.ReadFilterType;
import com.hartwig.hmftools.esvee.prep.types.ReadGroup;
import com.hartwig.hmftools.esvee.prep.types.ReadGroupStatus;
import com.hartwig.hmftools.esvee.common.ReadIdTrimmer;
import com.hartwig.hmftools.esvee.prep.types.PrepRead;
import com.hartwig.hmftools.esvee.prep.types.ReadType;
import com.hartwig.hmftools.esvee.prep.types.RemoteJunction;
import com.hartwig.hmftools.esvee.common.IndelCoords;

public class JunctionTracker
{
    private final ChrBaseRegion mRegion;
    private final PrepConfig mConfig;
    private final ReadFilterConfig mFilterConfig;
    private final List<BaseRegion> mBlacklistRegions;
    private final List<KnownHotspot> mKnownHotspots;
    private final BlacklistLocations mBlacklist;

    private final Map<String,ReadGroup> mReadGroupMap; // keyed by readId
    private final Set<String> mExpectedReadIds; // as indicated by another partition
    private final List<ReadGroup> mExpectedReadGroups;

    private final DiscordantGroups mDiscordantGroupFinder;

    // reads with their mate(s) in another partition, may or may not end up supporting a local junction
    private final Set<ReadGroup> mRemoteCandidateReadGroups;

    private final List<ReadGroup> mCandidateDiscordantGroups;

    private final List<JunctionData> mJunctions; // ordered by position
    private int mLastJunctionIndex;

    private ReadIdTrimmer mReadIdTrimmer;
    private int mInitialSupportingFrags;
    private final DiscordantStats mDiscordantStats;

    private final List<PerformanceCounter> mPerfCounters;

    private enum PerfCounters
    {
        InitJunctions,
        JunctionSupport,
        DiscordantGroups,
        JunctionFilter;
    };

    public JunctionTracker(
            final ChrBaseRegion region, final PrepConfig config, final HotspotCache hotspotCache, final BlacklistLocations blacklist)
    {
        mRegion = region;
        mConfig = config;
        mFilterConfig = config.ReadFiltering.config();

        mKnownHotspots = hotspotCache.findRegionHotspots(region);

        mDiscordantGroupFinder = new DiscordantGroups(mRegion, mFilterConfig.observedFragLengthMax(), mKnownHotspots, mConfig.TrackRemotes);

        mBlacklist = blacklist;
        mBlacklistRegions = Lists.newArrayList();

        List<BaseRegion> chrRegions = blacklist.getRegions(mRegion.Chromosome);

        if(chrRegions != null)
        {
            // extract the blacklist regions just for this partition
            chrRegions.stream().filter(x -> positionsOverlap(mRegion.start(), mRegion.end(), x.start(), x.end()))
                    .forEach(x -> mBlacklistRegions.add(x));
        }

        mReadGroupMap = Maps.newHashMap();
        mExpectedReadIds = Sets.newHashSet();
        mExpectedReadGroups = Lists.newArrayList();
        mRemoteCandidateReadGroups = Sets.newHashSet();
        mCandidateDiscordantGroups = Lists.newArrayList();
        mJunctions = Lists.newArrayList();
        mLastJunctionIndex = -1;
        mInitialSupportingFrags = 0;
        mDiscordantStats = new DiscordantStats();
        mReadIdTrimmer = new ReadIdTrimmer(mConfig.TrimReadId);

        mPerfCounters = Lists.newArrayList();

        for(PerfCounters pc : PerfCounters.values())
        {
            mPerfCounters.add(pc.ordinal(), new PerformanceCounter(pc.toString()));
        }
    }

    public void clear()
    {
        mReadGroupMap.clear();
        mExpectedReadIds.clear();
        mExpectedReadGroups.clear();
        mRemoteCandidateReadGroups.clear();
        mCandidateDiscordantGroups.clear();
        mJunctions.clear();
    }

    public List<JunctionData> junctions() { return mJunctions; }
    public List<PerformanceCounter> perfCounters() { return mPerfCounters; }
    public DiscordantStats discordantStats() { return mDiscordantStats; }

    public List<ReadGroup> formUniqueAssignedGroups()
    {
        // fragments can be added to more than one junction, so now collect up the unique ones
        List<ReadGroup> junctionGroups = Lists.newArrayList();
        Set<String> readIds = Sets.newHashSet();

        for(JunctionData junction : mJunctions)
        {
            ReadGroup.addUniqueReadGroups(readIds, junctionGroups, junction.junctionGroups());
            ReadGroup.addUniqueReadGroups(readIds, junctionGroups, junction.supportingGroups());
            ReadGroup.addUniqueReadGroups(readIds, junctionGroups, junction.exactSupportGroups());
        }

        // also gather expected remote reads and mark them as such
        for(ReadGroup readGroup : mExpectedReadGroups)
        {
            if(readIds.contains(readGroup.id()))
                continue;

            readGroup.setGroupState(ReadGroupStatus.EXPECTED);
            readGroup.reads().forEach(x -> x.setReadType(ReadType.EXPECTED));
            junctionGroups.add(readGroup);
        }

        return junctionGroups;
    }

    public List<ReadGroup> getRemoteCandidateReadGroups()
    {
        // gather groups with a read in another partition and not linked to a junction
        // to then pass to the combined cache
        return mRemoteCandidateReadGroups.stream()
                .filter(x -> !mExpectedReadGroups.contains(x))
                .filter(x -> x.noRegisteredJunctionPositions()).collect(Collectors.toList());
    }

    public int initialSupportingFrags() { return mInitialSupportingFrags; }

    public void processRead(final PrepRead read)
    {
        if(!read.hasMate() && !read.hasSuppAlignment())
        {
            // handle unpaired reads similarly to simple groups
            if(read.readType() == NO_SUPPORT || readInBlacklist(read))
                return;
        }

        String readId = mReadIdTrimmer.trim(read.id());

        ReadGroup readGroup = mReadGroupMap.get(readId);

        if(readGroup == null)
        {
            readGroup = new ReadGroup(read, readId);
            mReadGroupMap.put(readGroup.id(), readGroup);
            return;
        }

        readGroup.addRead(read);

        if(readGroup.isSimpleComplete()) // purge irrelevant groups
        {
            if(readGroup.allNoSupport())
            {
                mReadGroupMap.remove(readGroup.id());
            }
            else if(groupInBlacklist(readGroup))
            {
                mReadGroupMap.remove(readGroup.id());
            }
        }
    }

    public void setExpectedReads(final Set<String> expectedReads) { mExpectedReadIds.addAll(expectedReads); }

    public void assignFragments()
    {
        assignJunctionFragmentsAndSupport();

        filterJunctions();

        findDiscordantGroups();
    }

    public void assignJunctionFragmentsAndSupport()
    {
        List<ReadGroup> candidateSupportGroups = Lists.newArrayList();

        perfCounterStart(PerfCounters.InitJunctions);

        int duplicateGroups = markSupplementaryDuplicates(mReadGroupMap, mReadIdTrimmer);

        if(duplicateGroups > 100)
        {
            SV_LOGGER.debug("region({}) marked {} supplementary duplicates", mRegion, duplicateGroups);
        }

        // create junctions from read groups and then assignment supporting of various kinds
        // NOTE: the read groups are not ordered by position until the discordant group routine below
        for(ReadGroup readGroup : mReadGroupMap.values())
        {
            if(readGroup.groupStatus() == ReadGroupStatus.DUPLICATE)
                continue;

            if(mExpectedReadIds.remove(readGroup.id()))
            {
                readGroup.markHasRemoteJunctionReads();
                readGroup.setGroupState(ReadGroupStatus.EXPECTED);
                mExpectedReadGroups.add(readGroup);
            }

            // ignore any group with a short overlapping fragment, likely adapter
            if(readGroup.reads().stream().anyMatch(x -> ReadFilterType.isSet(x.filters(), INSERT_MAP_OVERLAP)))
                continue;

            if(readGroup.allNoSupport()) // ignore groups with only fully-filtered reads
                continue;

            // read groups can be assigned to more than one junction
            if(readGroup.hasReadType(ReadType.JUNCTION))
            {
                createJunction(readGroup);
                candidateSupportGroups.add(readGroup);
            }
            else if(readGroup.hasReadType(ReadType.CANDIDATE_SUPPORT))
            {
                candidateSupportGroups.add(readGroup);
            }
        }

        mJunctions.forEach(x -> x.setInitialRead(LOW_BASE_QUAL_THRESHOLD));

        perfCounterStop(PerfCounters.InitJunctions);

        mLastJunctionIndex = -1; // reset before supporting fragment assignment
        mInitialSupportingFrags = candidateSupportGroups.size();

        perfCounterStart(PerfCounters.JunctionSupport);

        // cannot early exit even if there are no junctions since could miss capture of any candidate support for remote junctions
        boolean hasJunctions = !mJunctions.isEmpty();

        // order by first read's start position to assist with efficient junction look-up using the last junction index
        if(hasJunctions)
            Collections.sort(candidateSupportGroups, new ReadGroup.ReadGroupComparator());

        Map<JunctionData,ReadType> supportedJunctions = Maps.newHashMap();
        for(ReadGroup readGroup : candidateSupportGroups)
        {
            supportedJunctions.clear();

            boolean hasBlacklistedRead = false;

            for(PrepRead read : readGroup.reads())
            {
                // supporting reads cannot fall in blacklist regions, despite allowing junctions in them
                if(readInBlacklist(read))
                {
                    hasBlacklistedRead = true;
                    continue;
                }
                else if(readMateInBlacklist(read))
                {
                    hasBlacklistedRead = true;
                }

                if(hasJunctions)
                {
                    // allow reads that match a junction to be considered for support (eg exact) for other junctions
                    if(read.readType() == ReadType.CANDIDATE_SUPPORT || read.readType() == ReadType.JUNCTION)
                        checkJunctionSupport(readGroup, read, supportedJunctions);
                }
                else if(hasBlacklistedRead)
                {
                    break;
                }
            }

            if(!supportedJunctions.isEmpty())
            {
                // check support from most important to least
                for(Map.Entry<JunctionData,ReadType> entry : supportedJunctions.entrySet())
                {
                    JunctionData junctionData = entry.getKey();

                    // group may already have a junction read, so skip if this is the case
                    if(readGroup.hasJunctionPosition(junctionData))
                        continue;

                    if(entry.getValue() == ReadType.EXACT_SUPPORT)
                        junctionData.addExactSupportGroup(readGroup);
                    else
                        junctionData.addSupportingGroup(readGroup);

                    readGroup.addJunctionPosition(junctionData);
                }
            }
            else if(!readGroup.hasJunctionPositions())
            {
                // will be passed onto the spanning cache if not assigned to a unfiltered junction
                mRemoteCandidateReadGroups.add(readGroup);
            }

            if(hasBlacklistedRead)
                continue;

            if(mDiscordantGroupFinder.isDiscordantGroup(readGroup))
            {
                mDiscordantStats.processReadGroup(readGroup);

                if(mDiscordantGroupFinder.isRelevantDiscordantGroup(readGroup))
                    mCandidateDiscordantGroups.add(readGroup);
            }
        }

        perfCounterStop(PerfCounters.JunctionSupport);
    }

    public void findDiscordantGroups()
    {
        if(mConfig.unpairedReads())
            return;

        perfCounterStart(PerfCounters.DiscordantGroups);

        List<JunctionData> discordantJunctions = mDiscordantGroupFinder.formDiscordantJunctions(mCandidateDiscordantGroups);

        if(mCandidateDiscordantGroups.size() > 2000 && !discordantJunctions.isEmpty())
        {
            SV_LOGGER.debug("region({}) found {} discordant group junctions from {} read groups",
                    mRegion, discordantJunctions.size(), mCandidateDiscordantGroups.size());
        }

        discordantJunctions.forEach(x -> addJunction(x));

        // no obvious need to re-check support at these junctions since all proximate facing read groups have already been tested
        // and allocated to these groups

        perfCounterStop(PerfCounters.DiscordantGroups);
    }

    private void createJunction(final ReadGroup readGroup)
    {
        List<JunctionData> junctions = Lists.newArrayListWithExpectedSize(2);
        List<RemoteJunction> remoteJunctions = mConfig.TrackRemotes ? Lists.newArrayList() : Collections.emptyList();

        for(PrepRead read : readGroup.reads())
        {
            if(read.isUnmapped())
                continue;

            if(mConfig.TrackRemotes && read.hasSuppAlignment())
            {
                RemoteJunction.addRemoteJunction(remoteJunctions, RemoteJunction.fromSupplementaryData(read.supplementaryAlignment()));
            }

            if(read.readType() != ReadType.JUNCTION)
                continue;

            IndelCoords indelCoords = read.indelCoords();

            if(indelCoords != null)
                handleIndelJunction(readGroup, read, indelCoords);

            if(ReadFilterType.isSet(read.filters(), SOFT_CLIP_LENGTH))
                continue;

            for(int i = 0; i <= 1; ++i)
            {
                int scLength = (i == 0) ? read.leftClipLength() : read.rightClipLength();

                // check with the shorter LINE soft-clip length since the soft-clip filter has already been checked, which takes LINE into account
                if(scLength < MIN_LINE_SOFT_CLIP_LENGTH)
                    continue;

                Orientation orientation = (i == 0) ? REVERSE : FORWARD;

                int position = orientation.isReverse() ? read.start() : read.end();

                // junctions cannot fall in blacklist regions
                if(positionInBlacklist(position))
                    continue;

                if(!mRegion.containsPosition(position))
                {
                    if(mConfig.TrackRemotes)
                        RemoteJunction.addRemoteJunction(remoteJunctions, new RemoteJunction(mRegion.Chromosome, position, orientation));
                }
                else
                {
                    JunctionData junctionData = getOrCreateJunction(read, orientation);
                    junctionData.addReadType(read, ReadType.JUNCTION);

                    if(!junctions.contains(junctionData))
                        junctions.add(junctionData);
                }
            }
        }

        if(junctions.isEmpty())
            return;

        junctions.forEach(x -> x.addJunctionReadGroup(readGroup));
        junctions.forEach(x -> readGroup.addJunctionPosition(x));

        for(RemoteJunction remoteJunction : remoteJunctions)
        {
            for(JunctionData junctionData : junctions)
            {
                // ignore remotes (typically) supplementaries which point to junction in another read in this group
                if(remoteJunction.matches(mRegion.Chromosome, junctionData.Position, junctionData.Orient))
                    continue;

                junctionData.addRemoteJunction(remoteJunction);
            }
        }
    }

    private boolean groupInBlacklist(final ReadGroup readGroup)
    {
        // test whether every read is in a blacklist region
        return readGroup.reads().stream().allMatch(x -> readInBlacklist(x));
    }

    private boolean readInBlacklist(final PrepRead read)
    {
        return mBlacklistRegions.stream().anyMatch(x -> positionsOverlap(x.start(), x.end(), read.start(), read.end()));
    }

    private boolean readMateInBlacklist(final PrepRead read)
    {
        if(!read.hasMate())
            return false;

        return mBlacklist.inBlacklistLocation(read.MateChromosome, read.MatePosStart, read.MatePosStart + mConfig.readLength());
    }

    private boolean positionInBlacklist(int junctionPosition)
    {
        return mBlacklistRegions.stream().anyMatch(x -> x.containsPosition(junctionPosition));
    }

    private void handleIndelJunction(final ReadGroup readGroup, final PrepRead read, final IndelCoords indelCoords)
    {
        if(indelCoords.Length < mFilterConfig.MinIndelLength)
            return;

        if(positionInBlacklist(indelCoords.PosStart) || positionInBlacklist(indelCoords.PosEnd))
            return;

        // a bit inefficient to search twice, but there won't be too many of these long indel reads
        JunctionData junctionStart = getOrCreateJunction(read, indelCoords.PosStart, FORWARD);
        JunctionData junctionEnd = getOrCreateJunction(read, indelCoords.PosEnd, REVERSE);

        junctionStart.markInternalIndel();
        junctionStart.addJunctionReadGroup(readGroup);
        junctionStart.addReadType(read, ReadType.JUNCTION);
        readGroup.addJunctionPosition(junctionStart);

        junctionEnd.markInternalIndel();
        junctionEnd.addJunctionReadGroup(readGroup);
        junctionEnd.addReadType(read, ReadType.JUNCTION);
        readGroup.addJunctionPosition(junctionEnd);
    }

    private void checkIndelSupport(final PrepRead read, final Map<JunctionData,ReadType> supportedJunctions)
    {
        IndelCoords indelCoords = read.indelCoords();

        if(indelCoords == null)
            return;

        int impliedUnclippedStart = read.start();
        int impliedUnclippedEnd = read.end();

        if(indelCoords.isInsert())
        {
            impliedUnclippedStart -= indelCoords.Length;
            impliedUnclippedEnd += indelCoords.Length;
        }
        else
        {
            impliedUnclippedStart += indelCoords.Length;
            impliedUnclippedEnd -= indelCoords.Length;
        }

        int readBoundsMin = min(read.start(), impliedUnclippedStart);
        int readBoundsMax = max(read.end(), impliedUnclippedEnd);

        // reads with a sufficiently long indel only need to cover a junction with any of their read bases, not the indel itself
        for(JunctionData junctionData : mJunctions)
        {
            if(min(readBoundsMin, readBoundsMax) > junctionData.Position)
                continue;

            if(junctionData.Position > max(readBoundsMin, readBoundsMax))
                break;

            if(supportedJunctions.containsKey(junctionData))
                continue;

            boolean hasExactSupport = false;

            for(int se = SE_START; se <= SE_END; ++se)
            {
                int indelPos = se == SE_START ? indelCoords.PosStart : indelCoords.PosEnd;

                if(indelPos != junctionData.Position)
                    continue;

                if(se == SE_START && junctionData.isReverse())
                    continue;
                else if(se == SE_END && junctionData.isForward())
                    continue;

                // indel coords support a junction
                read.setReadType(ReadType.EXACT_SUPPORT, true);
                junctionData.addReadType(read, ReadType.EXACT_SUPPORT);
                supportedJunctions.put(junctionData, ReadType.EXACT_SUPPORT);
                hasExactSupport = true;
            }

            if(!hasExactSupport)
            {
                read.setReadType(ReadType.SUPPORT, true);
                junctionData.addReadType(read, ReadType.SUPPORT);
                supportedJunctions.put(junctionData, ReadType.SUPPORT);
            }
        }
    }

    private JunctionData getOrCreateJunction(final PrepRead read, final Orientation orientation)
    {
        int junctionPosition = orientation.isReverse() ? read.start() : read.end();
        return getOrCreateJunction(read, junctionPosition, orientation);
    }

    private JunctionData getOrCreateJunction(final PrepRead read, final int junctionPosition, final Orientation orientation)
    {
        // junctions are stored in ascending order to make finding them more efficient, especially for supporting reads

        // first check the last index for a match
        if(mLastJunctionIndex >= 0)
        {
            JunctionData junctionData = mJunctions.get(mLastJunctionIndex);

            if(junctionData.Position == junctionPosition && junctionData.Orient == orientation)
                return junctionData;
        }

        int index = 0;

        while(index < mJunctions.size())
        {
            JunctionData junctionData = mJunctions.get(index);

            if(junctionData.Position == junctionPosition)
            {
                if(junctionData.Orient == orientation)
                {
                    setLastJunctionIndex(index);
                    return junctionData;
                }
            }
            else if(junctionData.Position > junctionPosition)
            {
                break;
            }

            ++index;
        }

        JunctionData junctionData = new JunctionData(junctionPosition, orientation, read);
        mJunctions.add(index, junctionData);
        setLastJunctionIndex(index);
        return junctionData;
    }

    private void addJunction(final JunctionData newJunction)
    {
        int index = 0;

        while(index < mJunctions.size())
        {
            JunctionData junctionData = mJunctions.get(index);

            if(junctionData.Position == newJunction.Position)
            {
                if(junctionData.Orient == newJunction.Orient)
                {
                    // favour discordant-only groups if the split-read junction has minimum support, to ensure if will be processed
                    // as a valid junction assembly by the assembler
                    if(!junctionData.discordantGroup() && newJunction.discordantGroup() && !junctionData.hotspot())
                    {
                        if(junctionData.junctionGroups().size() < mFilterConfig.MinJunctionSupport)
                        {
                            mJunctions.set(index, newJunction);
                        }
                    }

                    return;
                }
            }
            else if(junctionData.Position > newJunction.Position)
            {
                break;
            }

            ++index;
        }

        mJunctions.add(index, newJunction);
    }

    private void checkJunctionSupport(final ReadGroup readGroup, final PrepRead read, final Map<JunctionData,ReadType> supportedJunctions)
    {
        // first check indel support
        checkIndelSupport(read, supportedJunctions);

        int maxSupportDistance = mConfig.unpairedReads() ? UNPAIRED_READ_JUNCTION_DISTANCE : mFilterConfig.maxSupportingFragmentDistance();

        // first check the last index since the next read is likely to be close by
        int closeJunctionIndex = -1;

        if(mLastJunctionIndex >= 0)
        {
            JunctionData junctionData = mJunctions.get(mLastJunctionIndex);
            if(readWithinJunctionRange(read, junctionData, maxSupportDistance))
                closeJunctionIndex = mLastJunctionIndex;
        }

        if(closeJunctionIndex == -1)
            closeJunctionIndex = findJunctionIndex(read, maxSupportDistance);

        if(closeJunctionIndex < 0)
            return;

        setLastJunctionIndex(closeJunctionIndex);
        checkReadSupportsJunction(readGroup, read, mJunctions.get(closeJunctionIndex), supportedJunctions);

        // check up and down from this location
        for(int i = 0; i <= 1; ++i)
        {
            boolean searchUp = (i == 0);

            int index = searchUp ? closeJunctionIndex + 1 : closeJunctionIndex - 1;

            while(index >= 0 && index < mJunctions.size())
            {
                JunctionData junctionData = mJunctions.get(index);

                if(!readWithinJunctionRange(read, junctionData, maxSupportDistance))
                    break;

                checkReadSupportsJunction(readGroup, read, junctionData, supportedJunctions);

                if(searchUp)
                    ++index;
                else
                    --index;
            }
        }
    }

    private void checkReadSupportsJunction(
            final ReadGroup readGroup, final PrepRead read, final JunctionData junctionData,
            final Map<JunctionData,ReadType> supportedJunctions)
    {
        if(readGroup.hasJunctionPosition(junctionData))
            return;

        ReadType readType = supportedJunctions.get(junctionData);

        if(readType != ReadType.EXACT_SUPPORT && hasExactJunctionSupport(read, junctionData, mFilterConfig))
        {
            junctionData.addReadType(read, ReadType.EXACT_SUPPORT);
            read.setReadType(ReadType.EXACT_SUPPORT, true);
            supportedJunctions.put(junctionData, ReadType.EXACT_SUPPORT);
            return;
        }

        if(readType != ReadType.SUPPORT && !mConfig.unpairedReads() && hasOtherJunctionSupport(read, junctionData, mFilterConfig))
        {
            junctionData.addReadType(read, ReadType.SUPPORT);
            read.setReadType(ReadType.SUPPORT, true);
            supportedJunctions.put(junctionData, ReadType.SUPPORT);
        }
    }

    private void setLastJunctionIndex(int index)
    {
        mLastJunctionIndex = index;
    }

    private int findJunctionIndex(final PrepRead read, int maxSupportDistance)
    {
        if(mJunctions.size() <= 5)
        {
            for(int index = 0; index < mJunctions.size(); ++index)
            {
                if(readWithinJunctionRange(read, mJunctions.get(index), maxSupportDistance))
                    return index;
            }

            return -1;
        }

        // binary search on junctions if the collection gets too large
        int currentIndex = mJunctions.size() / 2;
        int lowerIndex = 0;
        int upperIndex = mJunctions.size() - 1;

        int iterations = 0;

        while(true)
        {
            JunctionData junctionData = mJunctions.get(currentIndex);

            if(readWithinJunctionRange(read, junctionData, maxSupportDistance))
                return currentIndex;

            if(read.end() < junctionData.Position)
            {
                // search lower
                if(lowerIndex + 1 == currentIndex)
                    return currentIndex;

                upperIndex = currentIndex;
                currentIndex = (lowerIndex + upperIndex) / 2;
            }
            else if(read.start() > junctionData.Position)
            {
                // search higher
                if(currentIndex + 1 == upperIndex)
                    return currentIndex;

                lowerIndex = currentIndex;
                currentIndex = (lowerIndex + upperIndex) / 2;
            }
            else
            {
                break;
            }

            ++iterations;

            if(iterations > 1000)
            {
                SV_LOGGER.warn("junction index search iterations({}}) junctions({}) index(cur={} low={} high={})",
                        iterations, mJunctions.size(), currentIndex, lowerIndex, upperIndex);
                break;
            }
        }

        return -1;
    }

    private boolean readWithinJunctionRange(final PrepRead read, final JunctionData junctionData, int maxDistance)
    {
        if(abs(read.end() - junctionData.Position) <= maxDistance)
            return true;

        if(abs(read.start() - junctionData.Position) <= maxDistance)
            return true;

        return false;
    }

    private void filterJunctions()
    {
        perfCounterStart(PerfCounters.JunctionFilter);

        Set<ReadGroup> removedReadGroups = Sets.newHashSet();

        // alternatively when processing then also just remove all processed junctions
        int index = 0;
        while(index < mJunctions.size())
        {
            JunctionData junctionData = mJunctions.get(index);

            if(junctionHasSupport(junctionData))
            {
                ++index;
            }
            else
            {
                mJunctions.remove(index);

                // now can remove candidate remote groups since they will be handled as part of actual junction groups
                junctionData.junctionGroups().forEach(x -> removedReadGroups.add(x));
                junctionData.exactSupportGroups().forEach(x -> removedReadGroups.add(x));
                junctionData.supportingGroups().forEach(x -> removedReadGroups.add(x));
            }
        }

        // reset read group junction positions, to remove those for purged junctions
        mReadGroupMap.values().forEach(x -> x.clearJunctionPositions());

        for(JunctionData junctionData : mJunctions)
        {
            junctionData.junctionGroups().forEach(x -> x.addJunctionPosition(junctionData));
            junctionData.exactSupportGroups().forEach(x -> x.addJunctionPosition(junctionData));
            junctionData.supportingGroups().forEach(x -> x.addJunctionPosition(junctionData));
        }

        // any reads no longer in any junction need to be reset to candidates only and will be passed to the spanning partition cache
        for(ReadGroup readGroup : removedReadGroups)
        {
            if(readGroup.groupStatus() == ReadGroupStatus.DUPLICATE)
            {
                readGroup.reads().forEach(x -> x.setReadType(NO_SUPPORT));
                continue;
            }

            if(!readGroup.hasJunctionPositions())
            {
                mRemoteCandidateReadGroups.add(readGroup);
                readGroup.reads().stream().filter(x -> x.readType() != NO_SUPPORT).forEach(x -> x.setReadType(ReadType.CANDIDATE_SUPPORT));
            }
        }

        perfCounterStop(PerfCounters.JunctionFilter);
    }

    private boolean junctionHasSupport(final JunctionData junctionData)
    {
        if(junctionData.discordantGroup())
            return true;

        // 1 junction read, 2 exact supporting reads altogether and 1 map-qual read
        int junctionFrags = junctionData.junctionFragmentCount();
        int exactSupportCount = junctionData.exactSupportFragmentCount();

        // check for a hotspot match
        if(junctionMatchesHotspot(mKnownHotspots, junctionData))
        {
            junctionData.markHotspot();

            if(junctionFrags + exactSupportCount >= MIN_HOTSPOT_JUNCTION_SUPPORT)
                return true;
        }

        if(!junctionData.internalIndel())
        {
            boolean hasPassingAlignedRead = false;

            for(PrepRead read : junctionData.readTypeReads().get(ReadType.JUNCTION))
            {
                if(aboveRepeatTrimmedAlignmentThreshold(read, mFilterConfig.MinCalcAlignmentScore, true))
                {
                    hasPassingAlignedRead = true;
                    break;
                }
            }

            if(!hasPassingAlignedRead)
                return false;
        }

        if(junctionFrags + exactSupportCount < mFilterConfig.MinJunctionSupport)
            return false;

        /*
        purgeSupplementaryDuplicates(junctionData);

        // re-assessed after duplicate assessment
        junctionFrags = junctionData.junctionFragmentCount();
        exactSupportCount = junctionData.exactSupportFragmentCount();
        */

        return junctionFrags + exactSupportCount >= mFilterConfig.MinJunctionSupport;
    }

    private void perfCounterStart(final PerfCounters pc)
    {
        if(!mConfig.PerfDebug)
            return;

        mPerfCounters.get(pc.ordinal()).start();
    }

    private void perfCounterStop(final PerfCounters pc)
    {
        if(!mConfig.PerfDebug)
            return;

        mPerfCounters.get(pc.ordinal()).stop();
    }
}
