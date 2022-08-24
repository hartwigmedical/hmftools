package com.hartwig.hmftools.svprep.reads;

import static java.lang.Math.abs;
import static java.lang.Math.max;
import static java.lang.Math.min;

import static com.hartwig.hmftools.common.utils.sv.BaseRegion.positionsOverlap;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_END;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_START;
import static com.hartwig.hmftools.common.utils.sv.SvCommonUtils.NEG_ORIENT;
import static com.hartwig.hmftools.common.utils.sv.SvCommonUtils.POS_ORIENT;
import static com.hartwig.hmftools.svprep.SvCommon.SV_LOGGER;
import static com.hartwig.hmftools.svprep.SvConstants.DEFAULT_READ_LENGTH;
import static com.hartwig.hmftools.svprep.SvConstants.LOW_BASE_QUALITY;
import static com.hartwig.hmftools.svprep.SvConstants.MAX_HIGH_QUAL_BASE_MISMATCHES;
import static com.hartwig.hmftools.svprep.SvConstants.MIN_HOTSPOT_JUNCTION_SUPPORT;
import static com.hartwig.hmftools.svprep.SvConstants.MIN_INDEL_SUPPORT_LENGTH;
import static com.hartwig.hmftools.svprep.SvConstants.MIN_LINE_SOFT_CLIP_LENGTH;
import static com.hartwig.hmftools.svprep.SvConstants.MIN_MAP_QUALITY;
import static com.hartwig.hmftools.svprep.reads.DiscordantGroups.formDiscordantJunctions;
import static com.hartwig.hmftools.svprep.reads.DiscordantGroups.isDiscordantGroup;
import static com.hartwig.hmftools.svprep.reads.ReadFilterType.INSERT_MAP_OVERLAP;
import static com.hartwig.hmftools.svprep.reads.ReadFilterType.POLY_G_SC;
import static com.hartwig.hmftools.svprep.reads.ReadFilterType.SOFT_CLIP_LENGTH;
import static com.hartwig.hmftools.svprep.reads.ReadFilters.isChimericRead;
import static com.hartwig.hmftools.svprep.reads.ReadGroup.addUniqueReadGroups;
import static com.hartwig.hmftools.svprep.reads.ReadRecord.findIndelCoords;
import static com.hartwig.hmftools.svprep.reads.ReadType.CANDIDATE_SUPPORT;
import static com.hartwig.hmftools.svprep.reads.ReadType.EXACT_SUPPORT;
import static com.hartwig.hmftools.svprep.reads.ReadType.EXPECTED;
import static com.hartwig.hmftools.svprep.reads.ReadType.JUNCTION;
import static com.hartwig.hmftools.svprep.reads.ReadType.SUPPORT;
import static com.hartwig.hmftools.svprep.reads.RemoteJunction.addRemoteJunction;

import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.stream.Collectors;

import com.beust.jcommander.internal.Sets;
import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.samtools.SoftClipSide;
import com.hartwig.hmftools.common.utils.PerformanceCounter;
import com.hartwig.hmftools.common.utils.sv.BaseRegion;
import com.hartwig.hmftools.common.utils.sv.ChrBaseRegion;
import com.hartwig.hmftools.svprep.BlacklistLocations;
import com.hartwig.hmftools.svprep.HotspotCache;
import com.hartwig.hmftools.svprep.SvConfig;

public class JunctionTracker
{
    private final ChrBaseRegion mRegion;
    private final SvConfig mConfig;
    private final ReadFilterConfig mFilterConfig;
    private final List<BaseRegion> mBlacklistRegions;
    private final List<ChrBaseRegion> mHotspotRegions;
    private final BlacklistLocations mBlacklist;

    private final Map<String,ReadGroup> mReadGroupMap; // keyed by readId
    private final Set<String> mExpectedReadIds; // as indicated by another partition
    private final List<ReadGroup> mExpectedReadGroups;
    private final List<ReadGroup> mRemoteCandidateReadGroups; // reads with their mate(s) in another partition, but not suppporting a junction
    private final List<ReadGroup> mCandidateDiscordantGroups;

    private final List<JunctionData> mJunctions; // ordered by position
    private int mLastJunctionIndex;

    private ReadIdTrimmer mReadIdTrimmer;
    private int mInitialSupportingFrags;
    private final int[] mBaseDepth;

    private final List<PerformanceCounter> mPerfCounters;

    private enum PerfCounters
    {
        InitJunctions,
        JunctionSupport,
        DiscordantGroups,
        JunctionFilter;
    };

    public JunctionTracker(
            final ChrBaseRegion region, final SvConfig config, final HotspotCache hotspotCache, final BlacklistLocations blacklist)
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
        mRemoteCandidateReadGroups = Lists.newArrayList();
        mCandidateDiscordantGroups = Lists.newArrayList();
        mJunctions = Lists.newArrayList();
        mLastJunctionIndex = -1;
        mInitialSupportingFrags = 0;
        mReadIdTrimmer = new ReadIdTrimmer(mConfig.TrimReadId);

        mBaseDepth = config.CaptureDepth ? new int[mRegion.baseLength()] : null;

        mPerfCounters = Lists.newArrayList();

        for(PerfCounters pc : PerfCounters.values())
        {
            mPerfCounters.add(pc.ordinal(), new PerformanceCounter(pc.toString()));
        }
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
            addUniqueReadGroups(readIds, junctionGroups, junction.JunctionGroups);
            addUniqueReadGroups(readIds, junctionGroups, junction.SupportingGroups);
            addUniqueReadGroups(readIds, junctionGroups, junction.ExactSupportGroups);
        }

        // also gather up expected remote reads and mark them as such
        for(ReadGroup readGroup : mExpectedReadGroups)
        {
            if(readIds.contains(readGroup.id()))
                continue;

            readGroup.setGroupState(ReadGroupStatus.EXPECTED);
            readGroup.reads().forEach(x -> x.setReadType(EXPECTED));
            junctionGroups.add(readGroup);
        }

        return junctionGroups;
    }

    public List<ReadGroup> getRemoteCandidateReadGroups()
    {
        // gather up groups with a read in another partition and not linked to a junction
        // to then pass to the combined cache
        return mRemoteCandidateReadGroups.stream().filter(x -> x.junctionPositions() == null).collect(Collectors.toList());
    }

    public int initialSupportingFrags() { return mInitialSupportingFrags; }

    public void processRead(final ReadRecord read)
    {
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
                captureDepth(readGroup);
                mReadGroupMap.remove(readGroup.id());
            }
            else if(groupInBlacklist(readGroup))
            {
                mReadGroupMap.remove(readGroup.id());
            }
        }
    }

    public void addExistingJunctions(final List<JunctionData> existingJunctions)
    {
        mJunctions.addAll(existingJunctions);
    }
    public void setExpectedReads(final Set<String> expectedReads) { mExpectedReadIds.addAll(expectedReads); }

    public void assignFragments()
    {
        assignJunctionFragmentsAndSupport();

        filterJunctions();

        findDiscordantGroups();

        if(mBaseDepth != null)
        {
            mJunctions.forEach(x -> x.setDepth(getBaseDepth(x.Position)));
        }
    }

    public void assignJunctionFragmentsAndSupport()
    {
        List<ReadGroup> candidateSupportGroups = Lists.newArrayList();

        perfCounterStart(PerfCounters.InitJunctions);

        // create junctions from read groups and then assignment supporting of various kinds
        // NOTE: the read groups are not ordered by position until the discordant group routine below
        for(ReadGroup readGroup : mReadGroupMap.values())
        {
            captureDepth(readGroup);

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
            if(readGroup.hasReadType(JUNCTION))
            {
                createJunction(readGroup);
                candidateSupportGroups.add(readGroup);
            }
            else if(readGroup.hasReadType(CANDIDATE_SUPPORT))
            {
                candidateSupportGroups.add(readGroup);
            }
        }

        mJunctions.forEach(x -> x.setInitialRead(mFilterConfig.MinSoftClipHighQual));

        perfCounterStop(PerfCounters.InitJunctions);

        if(mJunctions.isEmpty())
            return;

        mLastJunctionIndex = -1; // reset before supporting fragment assignment
        mInitialSupportingFrags = candidateSupportGroups.size();

        perfCounterStart(PerfCounters.JunctionSupport);

        // order by first read's start position to assist with efficient junction look-up using the last junction index
        Collections.sort(candidateSupportGroups, new ReadGroup.ReadGroupComparator());

        Map<JunctionData,ReadType> supportedJunctions = Maps.newHashMap();
        for(ReadGroup readGroup : candidateSupportGroups)
        {
            supportedJunctions.clear();

            boolean hasBlacklistedRead = false;

            for(ReadRecord read : readGroup.reads())
            {
                // supporting reads cannot fall in blacklist regions
                boolean readInBlacklist = readInBlacklist(read);

                if(readInBlacklist)
                {
                    hasBlacklistedRead = true;
                    // continue;
                }
                else if(readMateInBlacklist(read))
                {
                    hasBlacklistedRead = true;
                }

                // allow reads that match a junction to be considered for support (eg exact) for other junctions
                if(read.readType() == CANDIDATE_SUPPORT || read.readType() == JUNCTION)
                    checkJunctionSupport(readGroup, read, readInBlacklist, supportedJunctions);
            }

            if(!supportedJunctions.isEmpty())
            {
                for(Map.Entry<JunctionData,ReadType> entry : supportedJunctions.entrySet())
                {
                    JunctionData junctionData = entry.getKey();

                    // group may already have a junction read, so skip if this is the case
                    if(readGroup.hasJunctionPosition(junctionData.Position))
                        continue;

                    if(entry.getValue() == EXACT_SUPPORT)
                        junctionData.ExactSupportGroups.add(readGroup);
                    else
                        junctionData.SupportingGroups.add(readGroup);

                    readGroup.addJunctionPosition(junctionData.Position);
                }
            }
            else
            {
                mRemoteCandidateReadGroups.add(readGroup);
            }

            if(mConfig.FindDiscordantGroups && !hasBlacklistedRead &&
            isDiscordantGroup(readGroup, mFilterConfig.fragmentLengthMin(), mFilterConfig.fragmentLengthMax()))
            {
                mCandidateDiscordantGroups.add(readGroup);
            }
        }

        perfCounterStop(PerfCounters.JunctionSupport);
    }

    public void findDiscordantGroups()
    {
        if(!mConfig.FindDiscordantGroups)
            return;

        perfCounterStart(PerfCounters.DiscordantGroups);

        if(mCandidateDiscordantGroups.size() > 1000)
        {
            SV_LOGGER.debug("region({}) checking discordant groups from {} read groups", mRegion, mCandidateDiscordantGroups.size());
        }

        List<JunctionData> discordantJunctions = formDiscordantJunctions(mRegion, mCandidateDiscordantGroups, mFilterConfig.fragmentLengthMax());

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

        for(ReadRecord read : readGroup.reads())
        {
            if(read.isUnmapped())
                continue;

            if(mConfig.TrackRemotes && read.hasSuppAlignment())
            {
                addRemoteJunction(remoteJunctions, RemoteJunction.fromSupplementaryData(read.supplementaryAlignment()));
            }

            if(read.readType() != JUNCTION)
                continue;

            final int[] indelCoords = findIndelCoords(read, mFilterConfig.MinIndelLength);

            if(indelCoords != null)
                handleIndelJunction(readGroup, read, indelCoords);

            SoftClipSide scSide = SoftClipSide.fromCigar(read.cigar());

            if(scSide == null || ReadFilterType.isSet(read.filters(), SOFT_CLIP_LENGTH) || scSide.Length < MIN_LINE_SOFT_CLIP_LENGTH)
                continue;

            byte orientation = scSide.isLeft() ? NEG_ORIENT : POS_ORIENT;
            int position = scSide.isLeft() ? read.start() : read.end();

            // junctions cannot fall in blacklist regions
            if(positionInBlacklist(position))
                continue;

            if(!mRegion.containsPosition(position))
            {
                if(mConfig.TrackRemotes)
                    addRemoteJunction(remoteJunctions, new RemoteJunction(mRegion.Chromosome, position, orientation));
            }
            else
            {
                JunctionData junctionData = getOrCreateJunction(read, orientation);
                junctionData.addReadType(read, JUNCTION);

                if(!reachedFragmentCap(junctionData.junctionFragmentCount()) && !junctions.contains(junctionData))
                    junctions.add(junctionData);
            }
        }

        if(junctions.isEmpty())
            return;

        junctions.forEach(x -> x.JunctionGroups.add(readGroup));
        junctions.forEach(x -> readGroup.addJunctionPosition(x.Position));

        for(RemoteJunction remoteJunction : remoteJunctions)
        {
            for(JunctionData junctionData : junctions)
            {
                // ignore remotes (typically) supplementaries which point to junction in another read in this group
                if(remoteJunction.matches(mRegion.Chromosome, junctionData.Position, junctionData.Orientation))
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

    private boolean readInBlacklist(final ReadRecord read)
    {
        return mBlacklistRegions.stream().anyMatch(x -> positionsOverlap(x.start(), x.end(), read.start(), read.end()));
    }

    private boolean readMateInBlacklist(final ReadRecord read)
    {
        return mBlacklist.inBlacklistLocation(read.MateChromosome, read.MatePosStart, read.MatePosStart + DEFAULT_READ_LENGTH);
    }

    private boolean positionInBlacklist(int junctionPosition)
    {
        return mBlacklistRegions.stream().anyMatch(x -> x.containsPosition(junctionPosition));
    }

    private void handleIndelJunction(final ReadGroup readGroup, final ReadRecord read, final int[] indelCoords)
    {
        if(positionInBlacklist(indelCoords[SE_START]) || positionInBlacklist(indelCoords[SE_END]))
            return;

        // a bit inefficient to search twice, but there won't be too many of these long indel reads
        JunctionData junctionStart = getOrCreateJunction(read, indelCoords[SE_START], POS_ORIENT);
        JunctionData junctionEnd = getOrCreateJunction(read, indelCoords[SE_END], NEG_ORIENT);

        if(reachedFragmentCap(junctionStart.junctionFragmentCount()) || reachedFragmentCap(junctionEnd.junctionFragmentCount()))
            return;

        junctionStart.markInternalIndel();
        junctionStart.JunctionGroups.add(readGroup);
        junctionStart.addReadType(read, JUNCTION);
        readGroup.addJunctionPosition(indelCoords[SE_START]);

        junctionEnd.markInternalIndel();
        junctionEnd.JunctionGroups.add(readGroup);
        junctionStart.addReadType(read, JUNCTION);
        readGroup.addJunctionPosition(indelCoords[SE_END]);
    }

    private void checkIndelSupport(final ReadRecord read, final Map<JunctionData,ReadType> supportedJunctions)
    {
        final int[] indelCoords = findIndelCoords(read, MIN_INDEL_SUPPORT_LENGTH);

        if(indelCoords == null)
            return;

        for(JunctionData junctionData : mJunctions)
        {
            if(indelCoords[SE_START] > junctionData.Position)
                continue;

            if(junctionData.Position > indelCoords[SE_END])
                break;

            if(reachedFragmentCap(junctionData.supportingFragmentCount()))
                continue;

            if(supportedJunctions.containsKey(junctionData))
                continue;

            for(int se = SE_START; se <= SE_END; ++se)
            {
                int indelPos = indelCoords[se];

                if(indelPos != junctionData.Position)
                    continue;

                if(se == SE_START && junctionData.Orientation != POS_ORIENT)
                    continue;
                else if(se == SE_END && junctionData.Orientation != NEG_ORIENT)
                    continue;

                // indel coords support a junction
                read.setReadType(EXACT_SUPPORT, true);
                junctionData.addReadType(read, EXACT_SUPPORT);
                supportedJunctions.put(junctionData, EXACT_SUPPORT);
            }
        }
    }

    private JunctionData getOrCreateJunction(final ReadRecord read, final byte orientation)
    {
        int junctionPosition = orientation == NEG_ORIENT ? read.start() : read.end();
        return getOrCreateJunction(read, junctionPosition, orientation);
    }

    private JunctionData getOrCreateJunction(final ReadRecord read, final int junctionPosition, final byte orientation)
    {
        // junctions are stored in ascending order to make finding them more efficient, especially for supporting reads

        // first check the last index for a match
        if(mLastJunctionIndex >= 0)
        {
            JunctionData junctionData = mJunctions.get(mLastJunctionIndex);

            if(junctionData.Position == junctionPosition && junctionData.Orientation == orientation)
                return junctionData;
        }

        int index = 0;

        while(index < mJunctions.size())
        {
            JunctionData junctionData = mJunctions.get(index);

            if(junctionData.Position == junctionPosition)
            {
                if(junctionData.Orientation == orientation)
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
                if(junctionData.Orientation == newJunction.Orientation)
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

    private void checkJunctionSupport(
            final ReadGroup readGroup, final ReadRecord read, boolean readInBlacklist, final Map<JunctionData,ReadType> supportedJunctions)
    {
        // first check indel support
        checkIndelSupport(read, supportedJunctions);

        int maxSupportDistance = readInBlacklist ? mFilterConfig.MinSupportingReadDistance : mFilterConfig.maxSupportingFragmentDistance();

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
        checkReadSupportsJunction(readGroup, read, readInBlacklist, mJunctions.get(closeJunctionIndex), supportedJunctions);

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

                checkReadSupportsJunction(readGroup, read, readInBlacklist, junctionData, supportedJunctions);

                if(searchUp)
                    ++index;
                else
                    --index;
            }
        }
    }

    private void checkReadSupportsJunction(
            final ReadGroup readGroup, final ReadRecord read, boolean readInBlacklist, final JunctionData junctionData,
            final Map<JunctionData,ReadType> supportedJunctions)
    {
        if(reachedFragmentCap(junctionData.supportingFragmentCount())) // to limit processing
            return;

        if(readGroup.hasJunctionPosition(junctionData.Position))
            return;

        ReadType readType = supportedJunctions.get(junctionData);

        if(readType != EXACT_SUPPORT && hasExactJunctionSupport(read, junctionData, mFilterConfig))
        {
            junctionData.addReadType(read, EXACT_SUPPORT);
            read.setReadType(EXACT_SUPPORT, true);
            supportedJunctions.put(junctionData, EXACT_SUPPORT);
            return;
        }

        if(readType != SUPPORT && !readInBlacklist && hasDiscordantJunctionSupport(read, junctionData, mFilterConfig))
        {
            junctionData.addReadType(read, SUPPORT);
            read.setReadType(SUPPORT, true);
            supportedJunctions.put(junctionData, SUPPORT);
        }
    }

    private boolean reachedFragmentCap(final int fragments)
    {
        return mConfig.JunctionFragmentCap > 0 && fragments >= mConfig.JunctionFragmentCap;
    }

    private void setLastJunctionIndex(int index)
    {
        mLastJunctionIndex = index;
    }

    private int findJunctionIndex(final ReadRecord read, int maxSupportDistance)
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

    private boolean readWithinJunctionRange(final ReadRecord read, final JunctionData junctionData, int maxDistance)
    {
        if(abs(read.end() - junctionData.Position) <= maxDistance)
            return true;

        if(abs(read.start() - junctionData.Position) <= maxDistance)
            return true;

        return false;
    }

    public static boolean hasDiscordantJunctionSupport(
            final ReadRecord read, final JunctionData junctionData, final ReadFilterConfig filterConfig)
    {
        // correct orientation
        if(junctionData.Orientation != read.orientation())
            return false;

        // correct side of the junction
        int junctionDistance = 0;

        if(junctionData.Orientation == POS_ORIENT)
        {
            if(read.end() > junctionData.Position)
                return false;

            junctionDistance = abs(read.end() - junctionData.Position);
        }
        else
        {
            if(read.start() < junctionData.Position || abs(read.end() - junctionData.Position) > filterConfig.maxSupportingFragmentDistance())
                return false;

            junctionDistance = abs(read.start() - junctionData.Position);
        }

        // any soft-clipping on the correct side if close to the junction
        if(junctionDistance <= filterConfig.MinSupportingReadDistance)
        {
            if(junctionData.Orientation == POS_ORIENT && read.record().getCigar().isRightClipped())
                return true;

            if(junctionData.Orientation == NEG_ORIENT && read.record().getCigar().isLeftClipped())
                return true;
        }

        // otherwise can be distant if chimeric
        if(junctionDistance <= filterConfig.maxSupportingFragmentDistance())
            return isChimericRead(read.record(), filterConfig);

        return false;
    }

    public static boolean hasExactJunctionSupport(
            final ReadRecord read, final JunctionData junctionData, final ReadFilterConfig filterConfig)
    {
        if(!read.cigar().isLeftClipped() && !read.cigar().isRightClipped())
            return false;

        /*
        eg if there is a candidate read with 101M50S at base 1000, how does a supporting read need to align and match?
         In this case the soft clip is at 1050. Check all reads which have a soft clip on the same side within +/- 50 bases.
        eg. you may find a read with soft clip at 1047 with 10S.
        Check if the 10 bases match the last 3 M and first 7S of the original read allowing for low base qual mismatches.
         If they do then that counts as an additional read even though it is much shorter soft clip and does not exactly match the base.
         */
        final ReadRecord juncRead = junctionData.topJunctionRead();

        if(junctionData.Orientation == POS_ORIENT)
        {
            if(!read.cigar().isRightClipped())
                return false;

            int readRightPos = read.end();

            if(readRightPos == junctionData.Position)
                return true;

            if(juncRead == null)
                return false;

            // within 50 bases with exact sequence match in between the soft clip locations
            if(abs(readRightPos - junctionData.Position) > filterConfig.MinSupportingReadDistance)
                return false;

            // must also overlap the junction
            int scLength = read.cigar().getLastCigarElement().getLength();

            if(readRightPos + scLength < junctionData.Position)
                return false;

            int readLength = read.readBases().length();
            int readEndPosIndex = readLength - scLength - 1;

            int juncReadLength = juncRead.readBases().length();
            int juncReadScLength = juncRead.cigar().getLastCigarElement().getLength();
            int juncReadEndPosIndex = juncReadLength - juncReadScLength - 1;
            int endPosDiff = juncRead.end() - readRightPos;

            int junctionReadOffset = juncReadEndPosIndex - readEndPosIndex - endPosDiff;

            // test against all the read's right soft-clipped bases
            int highQualMismatches = 0;
            for(int i = readLength - scLength; i < readLength; ++i)
            {
                char readBase = read.readBases().charAt(i);

                int juncIndex = i + junctionReadOffset;
                if(juncIndex < 0 || juncIndex >= juncReadLength)
                    return false;

                char juncReadBase = juncRead.readBases().charAt(juncIndex);

                if(readBase == juncReadBase)
                    continue;

                if(read.baseQualities()[i] < LOW_BASE_QUALITY || juncRead.baseQualities()[juncIndex] < LOW_BASE_QUALITY)
                    continue;

                ++highQualMismatches;

                if(highQualMismatches > MAX_HIGH_QUAL_BASE_MISMATCHES)
                    return false;
            }
        }
        else
        {
            if(!read.cigar().isLeftClipped())
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

            int scLength = read.cigar().getFirstCigarElement().getLength();

            if(readLeftPos - scLength > junctionData.Position)
                return false;

            int juncReadScLength = juncRead.cigar().getFirstCigarElement().getLength();
            int posOffset = juncRead.start() - readLeftPos;
            int softClipDiff = juncReadScLength - scLength;
            int junctionReadOffset = softClipDiff - posOffset;
            int juncReadLength = juncRead.readBases().length();
            int highQualMismatches = 0;

            for(int i = 0; i < scLength; ++i)
            {
                char readBase = read.readBases().charAt(i);

                int juncIndex = i + junctionReadOffset;
                if(juncIndex < 0 || juncIndex >= juncReadLength)
                    return false;

                char juncReadBase = juncRead.readBases().charAt(juncIndex);

                if(readBase == juncReadBase)
                    continue;

                if(read.baseQualities()[i] < LOW_BASE_QUALITY || juncRead.baseQualities()[juncIndex] < LOW_BASE_QUALITY)
                    continue;

                ++highQualMismatches;

                if(highQualMismatches > MAX_HIGH_QUAL_BASE_MISMATCHES)
                    return false;
            }
        }

        return true;
    }

    private void filterJunctions()
    {
        perfCounterStart(PerfCounters.JunctionFilter);

        // alternatively when processing then also just remove all processed junctions
        int index = 0;
        while(index < mJunctions.size())
        {
            JunctionData junctionData = mJunctions.get(index);

            if(junctionHasSupport(junctionData))
                ++index;
            else
                mJunctions.remove(index);
        }

        // reset read group junction positions, to remove those for purged junctions
        mReadGroupMap.values().forEach(x -> x.clearJunctionPositions());

        for(JunctionData junctionData : mJunctions)
        {
            junctionData.JunctionGroups.forEach(x -> x.addJunctionPosition(junctionData.Position));
            junctionData.ExactSupportGroups.forEach(x -> x.addJunctionPosition(junctionData.Position));
            junctionData.SupportingGroups.forEach(x -> x.addJunctionPosition(junctionData.Position));
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

        boolean hasPassingMapQualRead = junctionData.ReadTypeReads.get(JUNCTION).stream().anyMatch(x -> x.mapQuality() > MIN_MAP_QUALITY);

        if(hasPassingMapQualRead && junctionFrags >= mFilterConfig.MinJunctionSupport)
            return true;

        // look in the exact matches for additional support
        if(!hasPassingMapQualRead)
        {
            hasPassingMapQualRead = junctionData.ReadTypeReads.get(EXACT_SUPPORT).stream().anyMatch(x -> x.mapQuality() > MIN_MAP_QUALITY);
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

    private void captureDepth(final ReadGroup readGroup)
    {
        if(mBaseDepth == null)
            return;

        int readsPosMin = readGroup.reads().stream().mapToInt(x -> x.start()).min().orElse(0);
        int readsPosMax = readGroup.reads().stream().mapToInt(x -> x.end()).max().orElse(0);
        int baseStart = max(readsPosMin - mRegion.start(), 0);
        int baseEnd = min(readsPosMax - mRegion.start(), mBaseDepth.length - 1);
        for(int i = baseStart; i <= baseEnd; ++i)
        {
            ++mBaseDepth[i];
        }
    }

    private int getBaseDepth(int position)
    {
        int baseIndex = position - mRegion.start();
        return baseIndex >= 0 && baseIndex < mBaseDepth.length ? mBaseDepth[baseIndex] : 0;
    }
}
