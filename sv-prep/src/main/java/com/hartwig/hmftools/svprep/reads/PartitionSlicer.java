package com.hartwig.hmftools.svprep.reads;

import static java.lang.String.format;

import static com.hartwig.hmftools.common.sv.ExcludedRegions.getPolyGRegion;
import static com.hartwig.hmftools.common.utils.PerformanceCounter.nanosToSeconds;
import static com.hartwig.hmftools.svprep.SvCommon.SV_LOGGER;
import static com.hartwig.hmftools.svprep.SvConstants.DOWN_SAMPLE_FRACTION;
import static com.hartwig.hmftools.svprep.SvConstants.DOWN_SAMPLE_THRESHOLD;
import static com.hartwig.hmftools.svprep.WriteType.BAM;
import static com.hartwig.hmftools.svprep.WriteType.READS;
import static com.hartwig.hmftools.svprep.reads.ReadType.CANDIDATE_SUPPORT;
import static com.hartwig.hmftools.svprep.reads.ReadType.EXPECTED;
import static com.hartwig.hmftools.svprep.reads.ReadType.JUNCTION;
import static com.hartwig.hmftools.svprep.reads.ReadType.RECOVERED;

import java.io.File;
import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.samtools.BamSlicer;
import com.hartwig.hmftools.common.utils.PerformanceCounter;
import com.hartwig.hmftools.common.utils.sv.ChrBaseRegion;
import com.hartwig.hmftools.svprep.CombinedReadGroups;
import com.hartwig.hmftools.svprep.CombinedStats;
import com.hartwig.hmftools.svprep.ExistingJunctionCache;
import com.hartwig.hmftools.svprep.ResultsWriter;
import com.hartwig.hmftools.svprep.SvConfig;
import com.hartwig.hmftools.svprep.WriteType;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;

public class PartitionSlicer
{
    private final int mId;
    private final SvConfig mConfig;
    private final ChrBaseRegion mRegion;
    private final CombinedReadGroups mCombinedReadGroups;
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
    private final PerformanceCounter[] mPerCounters;

    private static final int PC_SLICE = 0;
    private static final int PC_JUNCTIONS = 1;
    private static final int PC_MATE_SLICE = 2;
    private static final int PC_TOTAL = 3;

    public PartitionSlicer(
            final int id, final ChrBaseRegion region, final SvConfig config, final CombinedReadGroups combinedReadGroups,
            final ExistingJunctionCache existingJunctionCache, final ResultsWriter writer, final CombinedStats combinedStats)
    {
        mId = id;
        mConfig = config;
        mReadFilters = config.ReadFiltering;
        mCombinedReadGroups = combinedReadGroups;
        mWriter = writer;
        mRegion = region;
        mCombinedStats = combinedStats;

        mJunctionTracker = new JunctionTracker(mRegion, mConfig, mConfig.Hotspots, mConfig.Blacklist);

        mJunctionTracker.addExistingJunctions(existingJunctionCache.getRegionJunctions(mRegion));

        mSamReader = mConfig.BamFile != null ?
                SamReaderFactory.makeDefault().referenceSequence(new File(mConfig.RefGenomeFile)).open(new File(mConfig.BamFile)) : null;

        mBamSlicer = new BamSlicer(0, false, true, false);
        mBamSlicer.setKeepUnmapped();

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

        mPerCounters = new PerformanceCounter[PC_TOTAL+1];
        mPerCounters[PC_SLICE] = new PerformanceCounter("Slice");
        mPerCounters[PC_JUNCTIONS] = new PerformanceCounter("Junctions");
        mPerCounters[PC_MATE_SLICE] = new PerformanceCounter("UnmatchedSlice");
        mPerCounters[PC_TOTAL] = new PerformanceCounter("Total");

        mLogReadIds = !mConfig.LogReadIds.isEmpty();
    }

    public void run()
    {
        SV_LOGGER.debug("processing region({})", mRegion);

        mPerCounters[PC_TOTAL].start();

        mPerCounters[PC_SLICE].start();

        mBamSlicer.slice(mSamReader, Lists.newArrayList(mRegion), this::processSamRecord);

        mPerCounters[PC_SLICE].stop();

        mPerCounters[PC_JUNCTIONS].start();

        mJunctionTracker.setExpectedReads(mCombinedReadGroups.getExpectedReadIds(mRegion));

        mJunctionTracker.createJunctions();
        mJunctionTracker.filterJunctions();

        mPerCounters[PC_JUNCTIONS].stop();

        writeData();

        if(mStats.TotalReads > 0)
        {
            SV_LOGGER.debug("region({}) complete, stats({})", mRegion, mStats.toString());
            SV_LOGGER.debug("region({}) filters({})", mRegion, ReadFilterType.filterCountsToString(mStats.ReadFilterCounts));
        }

        mPerCounters[PC_TOTAL].stop();

        mCombinedStats.addPartitionStats(mStats);
        mCombinedStats.addPerfCounters(mPerCounters);

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
            if(mFilterRegion.containsPosition(readStart) || mFilterRegion.containsPosition(readStart + mConfig.ReadLength))
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

        if(filters != 0)
        {
            // allow low map quality through at this stage
            if(filters != ReadFilterType.MIN_MAP_QUAL.flag())
            {
                processFilteredRead(record, filters);
                return;
            }
        }

        ReadRecord read = ReadRecord.from(record);
        read.setFilters(filters);
        read.setReadType(JUNCTION);

        mJunctionTracker.processRead(read);
    }

    private void processFilteredRead(final SAMRecord record, final int filters)
    {
        // check criteria to keep an otherwise filtered, to see if it supports a non-filtered read or location
        // record filters by type
        for(ReadFilterType type : ReadFilterType.values())
        {
            if(type.isSet(filters))
                ++mStats.ReadFilterCounts[type.index()];
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

        if(captureCompleteGroups)
        {
            // read groups that span chromosomes or partitions need to be complete, so gather up their state to enable this
            Map<String,ReadGroup> spanningGroupsMap = Maps.newHashMap();

            mJunctionTracker.junctionGroups().forEach(x -> assignReadGroup(x, spanningGroupsMap));
            mJunctionTracker.supportingGroups().forEach(x -> assignReadGroup(x, spanningGroupsMap));

            List<ReadGroup> expectedGroups = mJunctionTracker.expectedGroups();
            expectedGroups.forEach(x -> assignReadGroup(x, spanningGroupsMap));
            expectedGroups.forEach(x -> x.reads().forEach(y -> y.setReadType(EXPECTED)));

            int spanningGroups = spanningGroupsMap.size();
            int totalGroups = mJunctionTracker.junctionGroups().size() + mJunctionTracker.supportingGroups().size() + expectedGroups.size();

            final Map<String, List<ExpectedRead>> missedReadsMap = Maps.newHashMap();
            mCombinedReadGroups.processSpanningReadGroups(mRegion, spanningGroupsMap, missedReadsMap);

            mPerCounters[PC_MATE_SLICE].start();
            List<ReadGroup> recoveredReadGroups = findMissedReads(missedReadsMap);
            mPerCounters[PC_MATE_SLICE].stop();

            if(totalGroups == 0 && missedReadsMap.isEmpty())
                return;

            int matchedReads = spanningGroupsMap.values().stream().mapToInt(x -> (int)x.reads().stream().filter(y -> y.written()).count()).sum();

            SV_LOGGER.debug("region({}) readGroups({}) complete(local={} spanning={}) expectedNotFound({}) matchedReads({})",
                    mRegion, totalGroups, totalGroups - spanningGroups, spanningGroups, expectedGroups.size(), matchedReads);

            if(mConfig.WriteTypes.contains(BAM))
            {
                mWriter.writeBamRecords(mJunctionTracker.junctionGroups());
                mWriter.writeBamRecords(mJunctionTracker.supportingGroups());
                mWriter.writeBamRecords(expectedGroups);
                mWriter.writeBamRecords(recoveredReadGroups);
            }

            if(mConfig.WriteTypes.contains(READS))
            {
                mWriter.writeReadGroup(mJunctionTracker.junctionGroups());
                mWriter.writeReadGroup(mJunctionTracker.supportingGroups());
                mWriter.writeReadGroup(expectedGroups);
                mWriter.writeReadGroup(recoveredReadGroups);
            }
        }

        if(mConfig.WriteTypes.contains(WriteType.JUNCTIONS))
        {
            mWriter.writeJunctionData(mRegion.Chromosome, mJunctionTracker.junctions());
        }

        mStats.JunctionCount += mJunctionTracker.junctions().size();
        mStats.JunctionFragmentCount += mJunctionTracker.junctionGroups().size();
        mStats.SupportingFragmentCount += mJunctionTracker.supportingGroups().size();
        mStats.InitialSupportingFragmentCount += mJunctionTracker.initialSupportingFrags();
    }

    private void assignReadGroup(final ReadGroup readGroup, Map<String,ReadGroup> partialGroupsMap)
    {
        readGroup.setPartitionCount(mRegion, mConfig.PartitionSize);

        if(readGroup.spansPartitions() && !readGroup.onlySupplementaries())
        {
            // only register remote reads if they are not supplementaries - instead just let these be written locally and not reconciled
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

    private static final int MAX_MISSED_READ_DEPTH = 100000;

    private List<ReadGroup> findMissedReads(final Map<String,List<ExpectedRead>> missedReadsMap)
    {
        List<ReadGroup> readGroups = Lists.newArrayList();

        if(missedReadsMap.isEmpty())
            return readGroups;

        int missedReadCount = missedReadsMap.values().stream().mapToInt(x -> x.size()).sum();
        SV_LOGGER.trace("region({}) searching for {} missed reads", mRegion, missedReadCount);

        // ignore reads in blacklist locations
        int blacklistCount = 0;

        // form 2 lists both order by missed read position
        Map<Integer,List<ExpectedRead>> positionReads = Maps.newHashMap();
        Map<Integer,List<String>> positionReadIds = Maps.newHashMap();

        for(Map.Entry<String, List<ExpectedRead>> entry : missedReadsMap.entrySet())
        {
            String readId = entry.getKey();

            for(ExpectedRead missedRead : entry.getValue())
            {
                if(mFilterRegion != null)
                {
                    if(mFilterRegion.containsPosition(missedRead.Position) || mFilterRegion.containsPosition(missedRead.Position + mConfig.ReadLength))
                        continue;
                }

                boolean inBlacklist = mConfig.Blacklist.inBlacklistLocation(
                        missedRead.Chromosome, missedRead.Position, missedRead.Position + mConfig.ReadLength);

                if(inBlacklist)
                {
                    ++blacklistCount;

                    if(!mConfig.RetrieveBlacklistMates)
                        continue;
                }

                List<ExpectedRead> posReads = positionReads.get(missedRead.Position);
                if(posReads == null)
                {
                    positionReads.put(missedRead.Position, Lists.newArrayList(missedRead));
                    positionReadIds.put(missedRead.Position, Lists.newArrayList(readId));
                }
                else
                {
                    posReads.add(missedRead);
                    positionReadIds.get(missedRead.Position).add(readId);
                }
            }
        }

        for(Map.Entry<Integer,List<ExpectedRead>> entry : positionReads.entrySet())
        {
            int position = entry.getKey();
            List<ExpectedRead> missedReads = entry.getValue();
            List<String> missedReadIds = positionReadIds.get(position);
            int posMissedReadCount = missedReadIds.size();

            long startTime = System.nanoTime();

            findMissedReads(readGroups, position, missedReads, missedReadIds);

            double sliceTime = nanosToSeconds(startTime, System.nanoTime());

            if(sliceTime > 2)
            {
                SV_LOGGER.debug("slice time({}) for {} missed reads, location({}:{})",
                        format("%.3f", sliceTime), posMissedReadCount, mRegion.Chromosome, position);
            }
        }

        readGroups.forEach(x -> x.setGroupState());

        if(missedReadsMap.size() != readGroups.size())
        {
            SV_LOGGER.debug("region({}) missed reads({}) recovered({}) in-blacklist({})",
                    mRegion, missedReadsMap.size(), readGroups.size(), blacklistCount);
        }

        return readGroups;
    }

    private void findMissedReads(
            List<ReadGroup> readGroups, int position, final List<ExpectedRead> missedReads, final List<String> missedReadIds)
    {
        final SAMRecordIterator iterator = mSamReader.queryAlignmentStart(mRegion.Chromosome, position);

        int readCount = 0;
        while(iterator.hasNext())
        {
            final SAMRecord record = iterator.next();
            ++readCount;

            if(readCount >= MAX_MISSED_READ_DEPTH)
                break;

            for(int i = 0; i < missedReads.size(); ++i)
            {
                if(missedReadIds.get(i).equals(record.getReadName()))
                {
                    ExpectedRead missedRead = missedReads.get(i);

                    if(missedRead.FirstInPair != record.getFirstOfPairFlag())
                        continue;

                    if(missedRead.IsSupplementary != record.getSupplementaryAlignmentFlag())
                        continue;

                    ReadRecord read = ReadRecord.from(record);
                    read.setReadType(RECOVERED);
                    readGroups.add(new ReadGroup(read));

                    missedReads.remove(i);
                    missedReadIds.remove(i);
                    break;
                }
            }

            if(missedReads.isEmpty())
                break;
        }

        iterator.close();
    }

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
