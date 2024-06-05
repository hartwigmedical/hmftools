package com.hartwig.hmftools.esvee.prep;

import static java.lang.Math.abs;
import static java.lang.Math.max;
import static java.lang.Math.min;

import static com.hartwig.hmftools.common.genome.region.Orientation.FORWARD;
import static com.hartwig.hmftools.common.genome.region.Orientation.REVERSE;
import static com.hartwig.hmftools.common.region.BaseRegion.positionWithin;
import static com.hartwig.hmftools.common.region.BaseRegion.positionsOverlap;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_END;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_START;
import static com.hartwig.hmftools.esvee.AssemblyConfig.SV_LOGGER;
import static com.hartwig.hmftools.esvee.common.SvConstants.LOW_BASE_QUAL_THRESHOLD;
import static com.hartwig.hmftools.esvee.prep.PrepConstants.MAX_HIGH_QUAL_BASE_MISMATCHES;
import static com.hartwig.hmftools.esvee.prep.PrepConstants.MIN_EXACT_BASE_PERC;
import static com.hartwig.hmftools.esvee.prep.PrepConstants.MIN_HOTSPOT_JUNCTION_SUPPORT;
import static com.hartwig.hmftools.esvee.prep.PrepConstants.MIN_LINE_SOFT_CLIP_LENGTH;
import static com.hartwig.hmftools.esvee.prep.PrepConstants.MIN_MAP_QUALITY;
import static com.hartwig.hmftools.esvee.prep.PrepConstants.UNPAIRED_READ_JUNCTION_DISTANCE;
import static com.hartwig.hmftools.esvee.prep.types.ReadFilterType.INSERT_MAP_OVERLAP;
import static com.hartwig.hmftools.esvee.prep.types.ReadFilterType.POLY_G_SC;
import static com.hartwig.hmftools.esvee.prep.types.ReadFilterType.SOFT_CLIP_LENGTH;
import static com.hartwig.hmftools.esvee.prep.types.ReadFilters.aboveRepeatTrimmedAlignmentThreshold;
import static com.hartwig.hmftools.esvee.prep.types.ReadFilters.isChimericRead;
import static com.hartwig.hmftools.esvee.prep.types.ReadType.NO_SUPPORT;

import static htsjdk.samtools.CigarOperator.M;
import static htsjdk.samtools.CigarOperator.S;

import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.stream.Collectors;

import com.beust.jcommander.internal.Sets;
import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.bam.ClippedSide;
import com.hartwig.hmftools.common.genome.region.Orientation;
import com.hartwig.hmftools.common.utils.PerformanceCounter;
import com.hartwig.hmftools.common.region.BaseRegion;
import com.hartwig.hmftools.common.region.ChrBaseRegion;
import com.hartwig.hmftools.esvee.prep.types.JunctionData;
import com.hartwig.hmftools.esvee.prep.types.JunctionsConfig;
import com.hartwig.hmftools.esvee.prep.types.ReadFilterConfig;
import com.hartwig.hmftools.esvee.prep.types.ReadFilterType;
import com.hartwig.hmftools.esvee.prep.types.ReadGroup;
import com.hartwig.hmftools.esvee.prep.types.ReadGroupStatus;
import com.hartwig.hmftools.esvee.common.ReadIdTrimmer;
import com.hartwig.hmftools.esvee.prep.types.PrepRead;
import com.hartwig.hmftools.esvee.prep.types.ReadType;
import com.hartwig.hmftools.esvee.prep.types.RemoteJunction;
import com.hartwig.hmftools.esvee.common.IndelCoords;

import htsjdk.samtools.CigarElement;

public class JunctionTracker
{
    private final ChrBaseRegion mRegion;
    private final JunctionsConfig mConfig;
    private final ReadFilterConfig mFilterConfig;
    private final List<BaseRegion> mBlacklistRegions;
    private final List<ChrBaseRegion> mHotspotRegions;
    private final BlacklistLocations mBlacklist;

    private final Map<String, ReadGroup> mReadGroupMap; // keyed by readId
    private final Set<String> mExpectedReadIds; // as indicated by another partition
    private final List<ReadGroup> mExpectedReadGroups;

    // reads with their mate(s) in another partition, may or may not end up supporting a local junction
    private final Set<ReadGroup> mRemoteCandidateReadGroups;

    private final List<ReadGroup> mCandidateDiscordantGroups;

    private final List<JunctionData> mJunctions; // ordered by position
    private int mLastJunctionIndex;

    private ReadIdTrimmer mReadIdTrimmer;
    private int mInitialSupportingFrags;

    private final List<PerformanceCounter> mPerfCounters;

    private enum PerfCounters
    {
        InitJunctions,
        JunctionSupport,
        DiscordantGroups,
        JunctionFilter;
    };

    public JunctionTracker(
            final ChrBaseRegion region, final PrepConfig svConfig, final HotspotCache hotspotCache, final BlacklistLocations blacklist)
    {
        this(region, JunctionsConfig.from(svConfig), hotspotCache, blacklist);
    }

    public JunctionTracker(
            final ChrBaseRegion region, final JunctionsConfig config, final HotspotCache hotspotCache, final BlacklistLocations blacklist)
    {
        mRegion = region;
        mConfig = config;
        mFilterConfig = config.ReadFiltering.config();

        mHotspotRegions = hotspotCache.findMatchingRegions(region);

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

    public List<ReadGroup> formUniqueAssignedGroups()
    {
        // fragments can be added to more than one junction, so now collect up the unique ones
        List<ReadGroup> junctionGroups = Lists.newArrayList();
        Set<String> readIds = Sets.newHashSet();

        for(JunctionData junction : mJunctions)
        {
            ReadGroup.addUniqueReadGroups(readIds, junctionGroups, junction.JunctionGroups);
            ReadGroup.addUniqueReadGroups(readIds, junctionGroups, junction.SupportingGroups);
            ReadGroup.addUniqueReadGroups(readIds, junctionGroups, junction.ExactSupportGroups);
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
        return mRemoteCandidateReadGroups.stream().filter(x -> x.noRegisteredJunctionPositions()).collect(Collectors.toList());
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

        // create junctions from read groups and then assignment supporting of various kinds
        // NOTE: the read groups are not ordered by position until the discordant group routine below
        for(ReadGroup readGroup : mReadGroupMap.values())
        {
            if(mExpectedReadIds.contains(readGroup.id()))
            {
                readGroup.markHasRemoteJunctionReads();
                mExpectedReadGroups.add(readGroup);
                mExpectedReadIds.remove(readGroup.id());
            }

            // ignore any group with a short overlapping fragment, likely adapter
            if(readGroup.reads().stream().anyMatch(x -> ReadFilterType.isSet(x.filters(), INSERT_MAP_OVERLAP)))
                continue;

            if(readGroup.allNoSupport()) // ignore groups with only fully-filtered reads
                continue;

            // ignore any group with a poly-G insert
            if(readGroup.reads().stream().anyMatch(x -> ReadFilterType.isSet(x.filters(), POLY_G_SC)))
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

        mJunctions.forEach(x -> x.setInitialRead(mFilterConfig.MinSoftClipHighQual));

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
                        junctionData.ExactSupportGroups.add(readGroup);
                    else
                        junctionData.SupportingGroups.add(readGroup);

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

            if(!mHotspotRegions.isEmpty()
            && DiscordantGroups.isDiscordantGroup(readGroup, mFilterConfig.fragmentLengthMin(), mFilterConfig.fragmentLengthMax()))
            {
                // require one end of this candidate group to be in a hotspot read
                boolean hasHotspotMatch = false;

                for(PrepRead read : readGroup.reads())
                {
                    if(mHotspotRegions.stream().anyMatch(x -> x.overlaps(read.Chromosome, read.start(), read.end())))
                    {
                        hasHotspotMatch = true;
                        break;
                    }
                }

                if(hasHotspotMatch)
                    mCandidateDiscordantGroups.add(readGroup);
            }
        }

        perfCounterStop(PerfCounters.JunctionSupport);
    }

    public void findDiscordantGroups()
    {
        if(mConfig.UnpairedReads)
            return;

        perfCounterStart(PerfCounters.DiscordantGroups);

        if(mCandidateDiscordantGroups.size() > 1000)
        {
            SV_LOGGER.debug("region({}) checking discordant groups from {} read groups", mRegion, mCandidateDiscordantGroups.size());
        }

        List<JunctionData> discordantJunctions = DiscordantGroups.formDiscordantJunctions(
                mRegion, mCandidateDiscordantGroups, mFilterConfig.fragmentLengthMax());

        if(!discordantJunctions.isEmpty())
        {
            SV_LOGGER.debug("region({}) found {} discordant group junctions", mRegion, discordantJunctions.size());
            discordantJunctions.forEach(x -> addJunction(x));
        }

        // no obvious need to re-check support at these junctions since all proximate facing read groups have already been tested
        // and allocated to these groups

        perfCounterStop(PerfCounters.DiscordantGroups);
    }

    private void createJunction(final ReadGroup readGroup)
    {
        List<JunctionData> junctions = Lists.newArrayList();
        List<RemoteJunction> remoteJunctions = Lists.newArrayList();

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

            ClippedSide scSide = ClippedSide.fromCigar(read.cigar(), false);

            if(scSide == null || ReadFilterType.isSet(read.filters(), SOFT_CLIP_LENGTH) || scSide.Length < MIN_LINE_SOFT_CLIP_LENGTH)
                continue;

            Orientation orientation = scSide.isLeft() ? REVERSE : FORWARD;
            int position = scSide.isLeft() ? read.start() : read.end();

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

        if(junctions.isEmpty())
            return;

        junctions.forEach(x -> x.JunctionGroups.add(readGroup));
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

        return mBlacklist.inBlacklistLocation(read.MateChromosome, read.MatePosStart, read.MatePosStart + mConfig.ReadLength);
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
        junctionStart.JunctionGroups.add(readGroup);
        junctionStart.addReadType(read, ReadType.JUNCTION);
        readGroup.addJunctionPosition(junctionStart);

        junctionEnd.markInternalIndel();
        junctionEnd.JunctionGroups.add(readGroup);
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
                    return;
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

        int maxSupportDistance = mConfig.UnpairedReads ? UNPAIRED_READ_JUNCTION_DISTANCE : mFilterConfig.maxSupportingFragmentDistance();

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

        if(readType != ReadType.SUPPORT && !mConfig.UnpairedReads && hasOtherJunctionSupport(read, junctionData, mFilterConfig))
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

                if(read.baseQualities()[i] < LOW_BASE_QUAL_THRESHOLD || juncRead.baseQualities()[juncIndex] < LOW_BASE_QUAL_THRESHOLD)
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

                if(read.baseQualities()[i] < LOW_BASE_QUAL_THRESHOLD || juncRead.baseQualities()[juncIndex] < LOW_BASE_QUAL_THRESHOLD)
                    continue;

                ++highQualMismatches;

                if(highQualMismatches > MAX_HIGH_QUAL_BASE_MISMATCHES)
                    return false;
            }

            double baseMatchPerc = baseMatches / (double)endIndex;
            return baseMatchPerc > MIN_EXACT_BASE_PERC;
        }
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
                junctionData.JunctionGroups.forEach(x -> removedReadGroups.add(x));
                junctionData.ExactSupportGroups.forEach(x -> removedReadGroups.add(x));
                junctionData.SupportingGroups.forEach(x -> removedReadGroups.add(x));
            }
        }

        // reset read group junction positions, to remove those for purged junctions
        mReadGroupMap.values().forEach(x -> x.clearJunctionPositions());

        for(JunctionData junctionData : mJunctions)
        {
            junctionData.JunctionGroups.forEach(x -> x.addJunctionPosition(junctionData));
            junctionData.ExactSupportGroups.forEach(x -> x.addJunctionPosition(junctionData));
            junctionData.SupportingGroups.forEach(x -> x.addJunctionPosition(junctionData));
        }

        // any reads no longer in any junction need to be reset to candidates only and will be passed to the spanning partition cache
        for(ReadGroup readGroup : removedReadGroups)
        {
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
        // first deal with junctions loaded from another sample - keep these if they've found any possible support
        if(junctionData.isExisting())
            return junctionData.totalFragmentCount() > 0;

        if(junctionData.discordantGroup())
            return true;

        // 1 junction read, 2 exact supporting reads altogether and 1 map-qual read
        int junctionFrags = junctionData.JunctionGroups.size();
        int exactSupportCount = junctionData.ExactSupportGroups.size();

        // check for a hotspot match
        if(mHotspotRegions.stream().anyMatch(x -> x.containsPosition(junctionData.Position)))
        {
            junctionData.markHotspot();

            if(junctionFrags + exactSupportCount >= MIN_HOTSPOT_JUNCTION_SUPPORT)
                return true;
        }

        boolean hasPassingMapQualRead = false;
        boolean hasPassingAlignedRead = false;

        for(PrepRead read : junctionData.ReadTypeReads.get(ReadType.JUNCTION))
        {
            hasPassingAlignedRead |= aboveRepeatTrimmedAlignmentThreshold(read, mFilterConfig.MinAlignmentBases);

            hasPassingMapQualRead |= read.mapQuality() >= MIN_MAP_QUALITY;
        }

        if(!hasPassingAlignedRead)
            return false;

        if(hasPassingMapQualRead && junctionFrags >= mFilterConfig.MinJunctionSupport)
            return true;

        // look in the exact matches for additional support
        if(!hasPassingMapQualRead)
        {
            hasPassingMapQualRead = junctionData.ReadTypeReads.get(ReadType.EXACT_SUPPORT).stream().anyMatch(x -> x.mapQuality() > MIN_MAP_QUALITY);
        }

        if(!hasPassingMapQualRead)
            return false;

        if(junctionFrags + exactSupportCount >= mFilterConfig.MinJunctionSupport)
            return true;

        return false;
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
