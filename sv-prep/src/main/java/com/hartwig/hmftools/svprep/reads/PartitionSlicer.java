package com.hartwig.hmftools.svprep.reads;

import static com.hartwig.hmftools.common.sv.ExcludedRegions.getPolyGRegion;
import static com.hartwig.hmftools.common.utils.sv.BaseRegion.positionsOverlap;
import static com.hartwig.hmftools.svprep.SvCommon.SV_LOGGER;
import static com.hartwig.hmftools.svprep.SvConstants.DOWN_SAMPLE_FRACTION;
import static com.hartwig.hmftools.svprep.SvConstants.DOWN_SAMPLE_THRESHOLD;
import static com.hartwig.hmftools.svprep.reads.ReadType.CANDIDATE_SUPPORT;
import static com.hartwig.hmftools.svprep.reads.ReadType.JUNCTION;

import java.util.List;
import java.util.Map;
import java.util.Set;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.samtools.BamSlicer;
import com.hartwig.hmftools.common.utils.PerformanceCounter;
import com.hartwig.hmftools.common.utils.sv.ChrBaseRegion;
import com.hartwig.hmftools.svprep.CombinedStats;
import com.hartwig.hmftools.svprep.ExistingJunctionCache;
import com.hartwig.hmftools.svprep.ResultsWriter;
import com.hartwig.hmftools.svprep.SpanningReadCache;
import com.hartwig.hmftools.svprep.SvConfig;
import com.hartwig.hmftools.svprep.WriteType;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReader;

public class PartitionSlicer
{
    private final int mId;
    private final SvConfig mConfig;
    private final ChrBaseRegion mRegion;
    private final SpanningReadCache mSpanningReadCache;
    private final ResultsWriter mWriter;
    private final ReadFilters mReadFilters;
    private final ChrBaseRegion mFilterRegion;

    private final SamReader mSamReader;
    private final BamSlicer mBamSlicer;

    private final JunctionTracker mJunctionTracker;

    private final PartitionStats mStats;
    private final CombinedStats mCombinedStats;

    private final ReadRateTracker mReadRateTracker;
    private boolean mRateLimitTriggered;
    private boolean mLogReadIds;
    private final List<PerformanceCounter> mPerfCounters;

    private enum PerfCounters
    {
        Slice,
        Junctions,
        Total;
    }

    public PartitionSlicer(
            final int id, final ChrBaseRegion region, final SvConfig config, final SamReader samReader, final BamSlicer bamSlicer,
            final SpanningReadCache spanningReadCache, final ExistingJunctionCache existingJunctionCache, final ResultsWriter writer,
            final CombinedStats combinedStats)
    {
        mId = id;
        mConfig = config;
        mReadFilters = config.ReadFiltering;
        mSpanningReadCache = spanningReadCache;
        mWriter = writer;
        mRegion = region;
        mCombinedStats = combinedStats;

        mJunctionTracker = new JunctionTracker(mRegion, mConfig, mConfig.Hotspots, mConfig.Blacklist);

        mJunctionTracker.addExistingJunctions(existingJunctionCache.getRegionJunctions(mRegion));

        mSamReader = samReader;
        mBamSlicer = bamSlicer;

        if(mConfig.ApplyDownsampling)
        {
            int rateSegmentLength = mConfig.PartitionSize / DOWN_SAMPLE_FRACTION;
            int downsampleThreshold = DOWN_SAMPLE_THRESHOLD / DOWN_SAMPLE_FRACTION;
            mReadRateTracker = new ReadRateTracker(rateSegmentLength, mRegion.start(), downsampleThreshold);
        }
        else
        {
            mReadRateTracker = null;
        }

        mRateLimitTriggered = false;

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

        mBamSlicer.slice(mSamReader, Lists.newArrayList(mRegion), this::processSamRecord);

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

        if(mConfig.PerfDebug)
            mPerfCounters.addAll(mJunctionTracker.perfCounters());

        mCombinedStats.addPerfCounters(mPerfCounters);

        if(mRateLimitTriggered)
            System.gc();
    }

    private static final boolean LOG_READ_ONLY = false;

    private void processSamRecord(final SAMRecord record)
    {
        int readStart = record.getAlignmentStart();

        if(!mRegion.containsPosition(readStart))
            return;

        ++mStats.TotalReads;

        if(mStats.TotalReads > 0 && (mStats.TotalReads % 1_000_000) == 0)
            System.gc();

        if(mFilterRegion != null)
        {
            if(positionsOverlap(readStart, readStart + mConfig.ReadLength, mFilterRegion.start(), mFilterRegion.end()))
                return;
        }

        if(mConfig.MaxPartitionReads > 0 && mStats.TotalReads >= mConfig.MaxPartitionReads)
        {
            SV_LOGGER.warn("region({}) readCount({}) exceeds maximum, stopping slice", mRegion, mStats.TotalReads);
            mBamSlicer.haltProcessing();
            return;
        }

        if(mLogReadIds) // debugging only
        {
            if(mConfig.LogReadIds.contains(record.getReadName()))
                SV_LOGGER.debug("specific readId({}) unmapped({})", record.getReadName(), record.getReadUnmappedFlag());
            else if(LOG_READ_ONLY)
                return;
        }

        if(readStart > 0 && !checkReadRateLimits(readStart))
            return;

        int filters = mReadFilters.checkFilters(record);

        if(filters == 0 || filters == ReadFilterType.MIN_MAP_QUAL.flag()) // allow reads only filtered by low map quality through
        {
            ReadRecord read = ReadRecord.from(record);
            read.setFilters(filters);
            read.setReadType(JUNCTION);

            mJunctionTracker.processRead(read);
        }
        else
        {
            processFilteredRead(record, filters);
        }
    }

    private void processFilteredRead(final SAMRecord record, final int filters)
    {
        // check criteria to keep an otherwise filtered, to see if it supports a non-filtered read or location
        if(mConfig.PerfDebug)
        {
            for(ReadFilterType type : ReadFilterType.values())
            {
                if(type.isSet(filters))
                    ++mStats.ReadFilterCounts[type.index()];
            }
        }

        // check for any evidence of support for an SV
        boolean isSupportCandidate = mReadFilters.isCandidateSupportingRead(record, filters);

        if(!isSupportCandidate && !mConfig.writeReads())
            return;

        ReadRecord read = ReadRecord.from(record);
        read.setFilters(filters);

        if(isSupportCandidate)
            read.setReadType(CANDIDATE_SUPPORT);

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

            mWriter.writeReadGroup(junctionGroups);
        }

        if(mConfig.WriteTypes.contains(WriteType.JUNCTIONS))
        {
            mWriter.writeJunctionData(mRegion.Chromosome, mJunctionTracker.junctions());
        }

        mStats.JunctionCount += mJunctionTracker.junctions().size();
        int junctionFragments = (int)junctionGroups.stream().filter(x -> x.hasReadType(JUNCTION)).count();
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

    private boolean checkReadRateLimits(int positionStart)
    {
        if(mReadRateTracker == null)
            return true;

        boolean wasLimited = mReadRateTracker.isRateLimited();
        int lastSegementReadCount = mReadRateTracker.readCount();

        boolean handleRead = mReadRateTracker.handleRead(positionStart);

        if(wasLimited != mReadRateTracker.isRateLimited())
        {
            if(mReadRateTracker.isRateLimited())
            {
                SV_LOGGER.info("region({}) rate limited with read count({}) at position({})",
                        mRegion, lastSegementReadCount, positionStart);
                mRateLimitTriggered = true;
            }
            else
            {
                SV_LOGGER.info("region({}) rate limit cleared at position({}), last read count({})",
                        mRegion, positionStart, lastSegementReadCount);
            }
        }

        return handleRead;
    }
}
