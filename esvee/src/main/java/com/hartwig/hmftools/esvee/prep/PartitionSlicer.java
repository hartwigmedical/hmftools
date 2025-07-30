package com.hartwig.hmftools.esvee.prep;

import static com.hartwig.hmftools.common.region.ExcludedRegions.getPolyGRegion;
import static com.hartwig.hmftools.common.region.BaseRegion.positionsOverlap;
import static com.hartwig.hmftools.esvee.assembly.AssemblyConfig.SV_LOGGER;
import static com.hartwig.hmftools.esvee.assembly.LineUtils.hasLineTail;
import static com.hartwig.hmftools.esvee.prep.PrepConstants.BAM_RECORD_SAMPLE_ID_TAG;
import static com.hartwig.hmftools.esvee.prep.PrepConstants.DEPTH_WINDOW_SIZE;
import static com.hartwig.hmftools.esvee.prep.types.WriteType.PREP_JUNCTION;

import java.util.List;
import java.util.Map;
import java.util.Set;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.bam.BamSlicer;
import com.hartwig.hmftools.common.perf.PerformanceCounter;
import com.hartwig.hmftools.common.region.BaseRegion;
import com.hartwig.hmftools.common.region.ChrBaseRegion;
import com.hartwig.hmftools.esvee.prep.types.CombinedStats;
import com.hartwig.hmftools.esvee.prep.types.PartitionStats;
import com.hartwig.hmftools.esvee.prep.types.ReadFilterType;
import com.hartwig.hmftools.esvee.prep.types.ReadGroup;
import com.hartwig.hmftools.esvee.prep.types.ReadGroupStatus;
import com.hartwig.hmftools.esvee.prep.types.PrepRead;
import com.hartwig.hmftools.esvee.prep.types.ReadType;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReader;

public class PartitionSlicer
{
    private final int mId;
    private final PrepConfig mConfig;
    private final ChrBaseRegion mRegion;
    private final SpanningReadCache mSpanningReadCache;
    private final ResultsWriter mWriter;
    private final ReadFilters mReadFilters;
    private final ChrBaseRegion mFilterRegion;

    private final List<SamReader> mSamReaders;
    private final BamSlicer mBamSlicer;
    private String mCurrentSampleId;

    private final JunctionTracker mJunctionTracker;
    private final DepthTracker mDepthTracker;

    private final PartitionStats mStats;
    private final CombinedStats mCombinedStats;

    private boolean mLogReadIds;
    private final List<PerformanceCounter> mPerfCounters;

    private enum PerfCounters
    {
        Slice,
        Junctions,
        Total;
    }

    public PartitionSlicer(
            final int id, final ChrBaseRegion region, final PrepConfig config, final List<SamReader> samReaders, final BamSlicer bamSlicer,
            final SpanningReadCache spanningReadCache, final ResultsWriter writer, final CombinedStats combinedStats)
    {
        mId = id;
        mConfig = config;
        mReadFilters = config.ReadFiltering;
        mSpanningReadCache = spanningReadCache;
        mWriter = writer;
        mRegion = region;
        mCombinedStats = combinedStats;

        mDepthTracker = new DepthTracker(new BaseRegion(mRegion.start(), mRegion.end()), DEPTH_WINDOW_SIZE);
        mJunctionTracker = new JunctionTracker(mRegion, mConfig, mDepthTracker, mConfig.Hotspots, mConfig.Blacklist);

        mSamReaders = samReaders;
        mBamSlicer = bamSlicer;

        ChrBaseRegion excludedRegion = getPolyGRegion(mConfig.RefGenVersion);
        mFilterRegion = region.overlaps(excludedRegion) ? excludedRegion : null;

        mStats = new PartitionStats();

        mPerfCounters = Lists.newArrayListWithExpectedSize(3);

        for(PerfCounters pc : PerfCounters.values())
        {
            mPerfCounters.add(pc.ordinal(), new PerformanceCounter(pc.toString()));
        }

        mLogReadIds = !mConfig.LogReadIds.isEmpty();
    }

    public void run()
    {
        SV_LOGGER.debug("processing region({})", mRegion);

        perfCounterStart(PerfCounters.Total);
        perfCounterStart(PerfCounters.Slice);

        for(int i = 0; i < mConfig.SampleIds.size(); ++i)
        {
            SamReader samReader = mSamReaders.get(i);
            mCurrentSampleId = mConfig.SampleIds.get(i);
            mBamSlicer.slice(samReader, mRegion, this::processSamRecord);
        }

        perfCounterStop(PerfCounters.Slice);

        perfCounterStart(PerfCounters.Junctions);

        Set<String> expectedJunctionReadIds = mSpanningReadCache.getExpectedReadIds(mRegion);
        mJunctionTracker.setExpectedReads(expectedJunctionReadIds);

        mJunctionTracker.assignFragments();

        perfCounterStop(PerfCounters.Junctions);

        writeData();

        if(mStats.TotalReads > 0)
        {
            SV_LOGGER.debug("region({}) complete, stats({})", mRegion, mStats.toString());
        }

        perfCounterStop(PerfCounters.Total);

        mCombinedStats.addPartitionStats(mStats);
        mCombinedStats.addDiscordantStats(mJunctionTracker.discordantStats());

        if(mConfig.PerfDebug)
            mPerfCounters.addAll(mJunctionTracker.perfCounters());

        mCombinedStats.addPerfCounters(mPerfCounters);
    }

    private void processSamRecord(final SAMRecord record)
    {
        int readStart = record.getAlignmentStart();

        if(!mRegion.containsPosition(readStart))
            return;

        ++mStats.TotalReads;

        if(mFilterRegion != null)
        {
            if(positionsOverlap(readStart, readStart + mConfig.readLength(), mFilterRegion.start(), mFilterRegion.end()))
                return;
        }

        if(mLogReadIds && mConfig.LogReadIds.contains(record.getReadName())) // debugging only
        {
            SV_LOGGER.debug("specific readId({}) unmapped({})", record.getReadName(), record.getReadUnmappedFlag());
        }

        boolean hasLineTail = hasLineTail(record);

        if(mReadFilters.ignoreRead(record, hasLineTail))
            return;

        record.setAttribute(BAM_RECORD_SAMPLE_ID_TAG, mCurrentSampleId);

        PrepRead read = new PrepRead(record);

        if(hasLineTail)
            read.markLineTail();

        mReadFilters.checkFilters(read);

        if(mConfig.PerfDebug)
        {
            for(ReadFilterType type : ReadFilterType.values())
            {
                if(type.isSet(read.filters()))
                    ++mStats.ReadFilterCounts[type.index()];
            }
        }

        if(read.unfiltered() || read.filters() == ReadFilterType.MIN_MAP_QUAL.flag())
        {
            // allow reads only filtered by low map quality through
            read.setReadType(ReadType.JUNCTION);

            mJunctionTracker.processRead(read);
        }
        else
        {
            processFilteredRead(read);
        }

        mDepthTracker.processRead(record);
    }

    private void processFilteredRead(final PrepRead read)
    {
        // check criteria to keep an otherwise filtered, to see if it supports a non-filtered read or location
        if(mConfig.PerfDebug)
        {
            for(ReadFilterType type : ReadFilterType.values())
            {
                if(type.isSet(read.filters()))
                    ++mStats.ReadFilterCounts[type.index()];
            }
        }

        // check for any evidence of support for an SV
        boolean isSupportCandidate = mReadFilters.isCandidateSupportingRead(read);

        if(!isSupportCandidate && !mConfig.writeReads())
            return;

        if(isSupportCandidate)
            read.setReadType(ReadType.CANDIDATE_SUPPORT);

        mJunctionTracker.processRead(read);
    }

    private void writeData()
    {
        boolean captureCompleteGroups = mConfig.writeReads();
        List<ReadGroup> junctionGroups = mJunctionTracker.formUniqueAssignedGroups();

        if(captureCompleteGroups)
        {
            // read groups that span chromosomes or partitions need to be completed, so gather up their state to enable this
            Map<String,ReadGroup> spanningGroupsMap = Maps.newHashMap();

            junctionGroups.forEach(x -> assignReadGroup(x, spanningGroupsMap));

            // inform the cache of incomplete groups with a possible remote junction or support read to avoid reslicing for these later
            List<ReadGroup> remoteCandidateGroups = mJunctionTracker.getRemoteCandidateReadGroups();
            remoteCandidateGroups.forEach(x -> assignReadGroup(x, spanningGroupsMap));

            int remoteCandidateCount = remoteCandidateGroups.size();
            int spanningGroupCount = spanningGroupsMap.size();
            int totalGroupCount = junctionGroups.size() + remoteCandidateCount;
            int expectedGroupCount = (int)junctionGroups.stream().filter(x -> x.groupStatus() == ReadGroupStatus.EXPECTED).count();

            mSpanningReadCache.processSpanningReadGroups(mRegion, spanningGroupsMap);

            if(totalGroupCount == 0)
                return;

            SV_LOGGER.debug("region({}) readGroups({} spanning={} candidate={}) expected({})",
                    mRegion, totalGroupCount, spanningGroupCount, remoteCandidateCount, expectedGroupCount);

            mWriter.writeReadGroups(junctionGroups);
        }

        if(mConfig.WriteTypes.contains(PREP_JUNCTION))
        {
            mWriter.writeJunctionData(mRegion.Chromosome, mJunctionTracker.junctions());
        }

        mStats.JunctionCount += mJunctionTracker.junctions().size();
        int junctionFragments = (int)junctionGroups.stream().filter(x -> x.hasReadType(ReadType.JUNCTION)).count();
        mStats.JunctionFragmentCount += junctionFragments;
        mStats.SupportingFragmentCount += junctionGroups.size() - junctionFragments;
        mStats.InitialSupportingFragmentCount += mJunctionTracker.initialSupportingFrags();
    }

    private void assignReadGroup(final ReadGroup readGroup, Map<String,ReadGroup> partialGroupsMap)
    {
        readGroup.setPartitionCount(mRegion, mConfig.PartitionSize);

        if(readGroup.spansPartitions() && readGroup.groupStatus() != ReadGroupStatus.EXPECTED)
        {
            ++mStats.SpanningGroups;

            readGroup.setGroupState();
            partialGroupsMap.put(readGroup.id(), readGroup);
        }
        else
        {
            readGroup.setGroupState();

            if(readGroup.isComplete())
                ++mStats.LocalCompleteGroups;
            else
                ++mStats.LocalIncompleteGroups;
        }
    }

    private void perfCounterStart(final PerfCounters pc)
    {
        if(mConfig.PerfDebug)
            mPerfCounters.get(pc.ordinal()).start(mRegion.toString());
        else
            mPerfCounters.get(pc.ordinal()).start();
    }

    private void perfCounterStop(final PerfCounters pc) { mPerfCounters.get(pc.ordinal()).stop(); }
}
