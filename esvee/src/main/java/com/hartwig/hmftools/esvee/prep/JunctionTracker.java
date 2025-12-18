package com.hartwig.hmftools.esvee.prep;

import static java.lang.Math.abs;
import static java.lang.Math.max;
import static java.lang.Math.min;

import static com.hartwig.hmftools.common.genome.region.Orientation.FORWARD;
import static com.hartwig.hmftools.common.genome.region.Orientation.REVERSE;
import static com.hartwig.hmftools.common.region.BaseRegion.positionWithin;
import static com.hartwig.hmftools.common.region.BaseRegion.positionsOverlap;
import static com.hartwig.hmftools.common.sequencing.SbxBamUtils.SBX_MAX_DUPLICATE_DISTANCE;
import static com.hartwig.hmftools.common.sv.StartEndIterator.SE_END;
import static com.hartwig.hmftools.common.sv.StartEndIterator.SE_START;
import static com.hartwig.hmftools.esvee.assembly.AssemblyConfig.SV_LOGGER;
import static com.hartwig.hmftools.esvee.common.SvConstants.hasPairedReads;
import static com.hartwig.hmftools.esvee.common.SvConstants.isSbx;
import static com.hartwig.hmftools.esvee.prep.JunctionUtils.INVALID_JUNC_INDEX;
import static com.hartwig.hmftools.esvee.prep.JunctionUtils.SIMPLE_SEARCH_COUNT;
import static com.hartwig.hmftools.esvee.prep.JunctionUtils.findJunctionMatchIndex;
import static com.hartwig.hmftools.esvee.prep.JunctionUtils.hasExactJunctionSupport;
import static com.hartwig.hmftools.esvee.prep.JunctionUtils.hasOtherJunctionSupport;
import static com.hartwig.hmftools.esvee.prep.JunctionUtils.hasWellAnchoredRead;
import static com.hartwig.hmftools.esvee.prep.JunctionUtils.markSupplementaryDuplicates;
import static com.hartwig.hmftools.esvee.prep.JunctionUtils.readWithinJunctionRange;
import static com.hartwig.hmftools.esvee.prep.KnownHotspot.junctionMatchesHotspot;
import static com.hartwig.hmftools.esvee.prep.PrepConstants.DEPTH_MIN_CHECK;
import static com.hartwig.hmftools.esvee.prep.PrepConstants.DEPTH_MIN_SUPPORT_RATIO_DISCORDANT;
import static com.hartwig.hmftools.esvee.prep.PrepConstants.DEPTH_MIN_SUPPORT_RATIO;
import static com.hartwig.hmftools.esvee.prep.PrepConstants.MIN_HOTSPOT_JUNCTION_SUPPORT;
import static com.hartwig.hmftools.esvee.prep.PrepConstants.MIN_LINE_SOFT_CLIP_LENGTH;
import static com.hartwig.hmftools.esvee.prep.PrepConstants.UNPAIRED_READ_JUNCTION_DISTANCE;
import static com.hartwig.hmftools.esvee.prep.types.DiscordantStats.isDiscordantUnpairedReadGroup;
import static com.hartwig.hmftools.esvee.prep.types.ReadFilterType.INSERT_MAP_OVERLAP;
import static com.hartwig.hmftools.esvee.prep.types.ReadFilterType.SOFT_CLIP_LENGTH;
import static com.hartwig.hmftools.esvee.prep.types.ReadType.NO_SUPPORT;

import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.stream.Collectors;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.genome.region.Orientation;
import com.hartwig.hmftools.common.perf.PerformanceCounter;
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
    private final List<KnownHotspot> mKnownHotspots;
    private final DepthTracker mDepthTracker;

    private final Map<String,ReadGroup> mReadGroupMap; // keyed by readId
    private final Set<String> mExpectedReadIds; // as indicated by another partition
    private final List<ReadGroup> mExpectedReadGroups;

    private final DiscordantGroups mDiscordantGroupFinder;

    // reads with their mate(s) in another partition, may or may not end up supporting a local junction
    private final Set<ReadGroup> mRemoteCandidateReadGroups;

    private final List<ReadGroup> mCandidateDiscordantGroups;

    private final List<JunctionData> mJunctions; // ordered by position
    private int mRecentJunctionIndex;

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
            final ChrBaseRegion region, final PrepConfig config, final DepthTracker depthTracker, final HotspotCache hotspotCache)
    {
        mRegion = region;
        mConfig = config;
        mFilterConfig = config.ReadFiltering.config();

        mDepthTracker = depthTracker;
        mKnownHotspots = hotspotCache.findRegionHotspots(region);

        mDiscordantGroupFinder = new DiscordantGroups(mRegion, mFilterConfig.observedFragLengthMax(), mKnownHotspots, mConfig.TrackRemotes);

        mReadGroupMap = Maps.newHashMap();
        mExpectedReadIds = Sets.newHashSet();
        mExpectedReadGroups = Lists.newArrayList();
        mRemoteCandidateReadGroups = Sets.newHashSet();
        mCandidateDiscordantGroups = Lists.newArrayList();
        mJunctions = Lists.newArrayList();
        mRecentJunctionIndex = INVALID_JUNC_INDEX;
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
            if(read.readType() == NO_SUPPORT)
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

        int permittedPositionDiff = isSbx() ? SBX_MAX_DUPLICATE_DISTANCE : 0; // matching Redux

        int duplicateGroups = markSupplementaryDuplicates(mReadGroupMap, mReadIdTrimmer, permittedPositionDiff);

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

        mJunctions.forEach(x -> x.setInitialRead());

        perfCounterStop(PerfCounters.InitJunctions);

        mRecentJunctionIndex = INVALID_JUNC_INDEX; // reset before supporting fragment assignment
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

            for(PrepRead read : readGroup.reads())
            {
                if(hasJunctions)
                {
                    // allow reads that match a junction to be considered for support (eg exact) for other junctions
                    if(read.readType() == ReadType.CANDIDATE_SUPPORT || read.readType() == ReadType.JUNCTION)
                        checkJunctionSupport(readGroup, read, supportedJunctions);
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

            boolean isDiscordantGroup = false;

            if(hasPairedReads())
            {
                isDiscordantGroup = mDiscordantGroupFinder.isDiscordantGroup(readGroup);

                if(isDiscordantGroup && mDiscordantGroupFinder.isRelevantDiscordantGroup(readGroup))
                    mCandidateDiscordantGroups.add(readGroup);
            }

            if(isDiscordantGroup || isDiscordantUnpairedReadGroup(readGroup))
                mDiscordantStats.processReadGroup(readGroup);
        }

        perfCounterStop(PerfCounters.JunctionSupport);
    }

    public void findDiscordantGroups()
    {
        if(!hasPairedReads())
            return;

        perfCounterStart(PerfCounters.DiscordantGroups);

        List<JunctionData> discordantJunctions = mDiscordantGroupFinder.formDiscordantJunctions(mCandidateDiscordantGroups);

        if(mCandidateDiscordantGroups.size() > 2000 && !discordantJunctions.isEmpty())
        {
            SV_LOGGER.debug("region({}) found {} discordant group junctions from {} read groups",
                    mRegion, discordantJunctions.size(), mCandidateDiscordantGroups.size());
        }

        discordantJunctions.forEach(x -> addDiscordantJunction(x));

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

                int position = orientation.isReverse() ? read.AlignmentStart : read.AlignmentEnd;

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

    private void handleIndelJunction(final ReadGroup readGroup, final PrepRead read, final IndelCoords indelCoords)
    {
        if(indelCoords.Length < mFilterConfig.MinIndelLength)
            return;

        // a bit inefficient to search twice, but there won't be too many of these long indel reads
        JunctionData junctionStart = getOrCreateJunction(read, indelCoords.PosStart, FORWARD);
        JunctionData junctionEnd = getOrCreateJunction(read, indelCoords.PosEnd, REVERSE);

        junctionStart.setLinkedIndel(junctionEnd);
        junctionEnd.setLinkedIndel(junctionStart);

        junctionStart.markInternalIndel();
        junctionStart.addJunctionReadGroup(readGroup);
        junctionStart.addReadType(read, ReadType.JUNCTION);
        readGroup.addJunctionPosition(junctionStart);

        junctionEnd.markInternalIndel();
        junctionEnd.addJunctionReadGroup(readGroup);
        junctionEnd.addReadType(read, ReadType.JUNCTION);
        readGroup.addJunctionPosition(junctionEnd);
    }

    private void checkIndelSupport(final PrepRead read, final Map<JunctionData,ReadType> supportedJunctions, final int closeJunctionIndex)
    {
        IndelCoords indelCoords = read.indelCoords();

        if(indelCoords == null)
            return;

        int impliedUnclippedStart = read.AlignmentStart;
        int impliedUnclippedEnd = read.AlignmentEnd;

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

        int readBoundsMin = min(read.AlignmentStart, impliedUnclippedStart);
        int readBoundsMax = max(read.AlignmentEnd, impliedUnclippedEnd);

        // reads with a sufficiently long indel only need to cover a junction with any of their read bases, not the indel itself
        for(int i = 0; i <= 1; ++i)
        {
            boolean searchUp = (i == 0);

            int index = searchUp ? closeJunctionIndex : closeJunctionIndex - 1;

            while(index >= 0 && index < mJunctions.size())
            {
                JunctionData junctionData = mJunctions.get(index);

                if(searchUp && junctionData.Position > readBoundsMax)
                    break;

                if(!searchUp && readBoundsMin > junctionData.Position)
                    break;

                if(supportedJunctions.containsKey(junctionData))
                {
                    index += searchUp ? 1 : -1;
                    continue;
                }

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

                index += searchUp ? 1 : -1;
            }
        }
    }

    private JunctionData getOrCreateJunction(final PrepRead read, final Orientation orientation)
    {
        int junctionPosition = orientation.isReverse() ? read.AlignmentStart : read.AlignmentEnd;
        return getOrCreateJunction(read, junctionPosition, orientation);
    }

    private JunctionData getOrCreateJunction(final PrepRead read, final int junctionPosition, final Orientation orientation)
    {
        // junctions are stored in ascending order

        // first check the last used index for a match
        if(mRecentJunctionIndex >= 0)
        {
            JunctionData junctionData = mJunctions.get(mRecentJunctionIndex);

            if(junctionData.Position == junctionPosition && junctionData.Orient == orientation)
                return junctionData;
        }

        int closeJuncIndex = findJunctionIndex(junctionPosition); // returns a match or the preceding index

        if(closeJuncIndex != INVALID_JUNC_INDEX)
        {
            int exactJuncIndex = findJunctionMatchIndex(mJunctions, junctionPosition, orientation, closeJuncIndex);

            if(exactJuncIndex != INVALID_JUNC_INDEX)
            {
                setLastJunctionIndex(exactJuncIndex);
                return mJunctions.get(exactJuncIndex);
            }
        }

        int newJuncIndex = closeJuncIndex + 1;

        JunctionData junctionData = new JunctionData(junctionPosition, orientation, read);

        if(newJuncIndex >= 0 && newJuncIndex < mJunctions.size())
            mJunctions.add(newJuncIndex, junctionData);
        else
            mJunctions.add(junctionData);

        setLastJunctionIndex(newJuncIndex);
        return junctionData;
    }

    @VisibleForTesting
    public void addDiscordantJunction(final JunctionData discJunction)
    {
        int closeJuncIndex = findJunctionIndex(discJunction.Position);

        if(closeJuncIndex != INVALID_JUNC_INDEX)
        {
            int exactJuncIndex = findJunctionMatchIndex(mJunctions, discJunction.Position, discJunction.Orient, closeJuncIndex);

            if(exactJuncIndex != INVALID_JUNC_INDEX)
            {
                JunctionData junctionData = mJunctions.get(exactJuncIndex);

                if(!junctionData.discordantGroup() && discJunction.discordantGroup() && !junctionData.hotspot()
                && junctionData.junctionGroups().size() < mFilterConfig.MinJunctionSupport)
                {
                    mJunctions.set(exactJuncIndex, discJunction);
                }

                return;
            }
        }

        int newJuncIndex = closeJuncIndex + 1;
        mJunctions.add(newJuncIndex, discJunction);
    }

    private void checkJunctionSupport(final ReadGroup readGroup, final PrepRead read, final Map<JunctionData,ReadType> supportedJunctions)
    {
        boolean hasPairedReads = hasPairedReads();
        int distanceBuffer = hasPairedReads ? mFilterConfig.maxSupportingFragmentDistance() : UNPAIRED_READ_JUNCTION_DISTANCE;

        // first check the last index since the next read is likely to be close by
        int closeJuncIndex = INVALID_JUNC_INDEX;

        if(mRecentJunctionIndex >= 0)
        {
            JunctionData junctionData = mJunctions.get(mRecentJunctionIndex);

            if(readWithinJunctionRange(read, junctionData, distanceBuffer)
            || positionWithin(junctionData.Position, read.AlignmentStart, read.AlignmentEnd))
            {
                closeJuncIndex = mRecentJunctionIndex;
            }
            else if(mRecentJunctionIndex < mJunctions.size() - 1)
            {
                // the preceding junction will typically have been found from the read start, so check it if is still good to use
                JunctionData nextJunctionData = mJunctions.get(mRecentJunctionIndex + 1);

                if(positionsOverlap(read.AlignmentStart, read.AlignmentEnd, junctionData.Position, nextJunctionData.Position))
                    closeJuncIndex = mRecentJunctionIndex;
            }
        }

        if(closeJuncIndex == INVALID_JUNC_INDEX)
            closeJuncIndex = findJunctionIndex(read.AlignmentStart);

        if(closeJuncIndex == INVALID_JUNC_INDEX)
            closeJuncIndex = 0;

        setLastJunctionIndex(closeJuncIndex);

        // check up and down from this location
        for(int i = 0; i <= 1; ++i)
        {
            boolean searchUp = (i == 0);

            int index = searchUp ? closeJuncIndex : closeJuncIndex - 1;

            while(index >= 0 && index < mJunctions.size())
            {
                JunctionData junctionData = mJunctions.get(index);

                if(!hasPairedReads)
                {
                    if(searchUp && junctionData.Position >= read.UnclippedEnd)
                        break;

                    if(!searchUp && junctionData.Position < read.UnclippedStart)
                        break;

                    if(!readWithinJunctionRange(read, junctionData, distanceBuffer))
                    {
                        index += searchUp ? 1 : -1;
                        continue;
                    }
                }
                else
                {
                    if(!readWithinJunctionRange(read, junctionData, distanceBuffer))
                    {
                        // given the find routine returns a match or the preceding, must search at least 1 higher
                        if(!searchUp || index > closeJuncIndex)
                            break;
                    }
                }

                checkReadSupportsJunction(readGroup, read, junctionData, supportedJunctions);

                index += searchUp ? 1 : -1;
            }
        }

        checkIndelSupport(read, supportedJunctions, closeJuncIndex);
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

        if(readType != ReadType.SUPPORT && hasPairedReads() && hasOtherJunctionSupport(read, junctionData, mFilterConfig))
        {
            junctionData.addReadType(read, ReadType.SUPPORT);
            read.setReadType(ReadType.SUPPORT, true);
            supportedJunctions.put(junctionData, ReadType.SUPPORT);
        }
    }

    private void setLastJunctionIndex(int index)
    {
        mRecentJunctionIndex = index;
    }

    private int findJunctionIndex(int position)
    {
        return JunctionUtils.findJunctionIndex(mJunctions, position, SIMPLE_SEARCH_COUNT);
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

            junctionData.setDepth(mDepthTracker.calcDepth(junctionData.Position));

            boolean junctionHasSupport = junctionHasSupport(junctionData);

            if(!junctionHasSupport && junctionData.hasLinkedIndel())
                junctionHasSupport = junctionHasSupport(junctionData.linkedIndel());

            if(junctionHasSupport)
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
        if(!junctionData.hotspot() && !junctionAboveMinAF(junctionData))
            return false;

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
            if(!hasWellAnchoredRead(junctionData, mFilterConfig))
                return false;
        }

        if(junctionFrags + exactSupportCount < mFilterConfig.MinJunctionSupport)
            return false;

        return junctionFrags + exactSupportCount >= mFilterConfig.MinJunctionSupport;
    }

    private static final int DEPTH_MIN_READ_COUNT = 200;

    private boolean junctionAboveMinAF(final JunctionData junctionData)
    {
        // sites with an AF below 0.5% (1% for discordant) are filtered
        int regionDepth = junctionData.depth();

        if(regionDepth < DEPTH_MIN_CHECK)
            return true;

        double requiredSupportRatio = junctionData.discordantGroup() ? DEPTH_MIN_SUPPORT_RATIO_DISCORDANT : DEPTH_MIN_SUPPORT_RATIO;

        int junctionSupport = junctionData.junctionFragmentCount() + junctionData.exactSupportFragmentCount();

        if(junctionData.discordantGroup())
            junctionSupport += junctionData.supportingFragmentCount();

        double requiredSupport = regionDepth * requiredSupportRatio;

        return junctionSupport >= requiredSupport;
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

    @VisibleForTesting
    public JunctionData addJunctionData(final PrepRead read)
    {
        Orientation orientation = read.isLeftClipped() ? REVERSE : FORWARD;
        return getOrCreateJunction(read, orientation);
    }
}
