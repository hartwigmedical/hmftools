package com.hartwig.hmftools.svprep.reads;

import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion.V37;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion.V38;
import static com.hartwig.hmftools.svprep.SvCommon.SV_LOGGER;
import static com.hartwig.hmftools.svprep.SvConstants.DOWN_SAMPLE_FRACTION;
import static com.hartwig.hmftools.svprep.SvConstants.DOWN_SAMPLE_THRESHOLD;
import static com.hartwig.hmftools.svprep.SvConstants.EXCLUDED_REGION_1_REF_37;
import static com.hartwig.hmftools.svprep.SvConstants.EXCLUDED_REGION_1_REF_38;
import static com.hartwig.hmftools.svprep.WriteType.BAM;
import static com.hartwig.hmftools.svprep.WriteType.READS;
import static com.hartwig.hmftools.svprep.reads.ReadGroup.MAX_GROUP_READ_COUNT;
import static com.hartwig.hmftools.svprep.reads.ReadType.CANDIDATE_SUPPORT;
import static com.hartwig.hmftools.svprep.reads.ReadType.JUNCTION;

import java.io.File;
import java.util.Map;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.samtools.BamSlicer;
import com.hartwig.hmftools.common.utils.PerformanceCounter;
import com.hartwig.hmftools.common.utils.sv.ChrBaseRegion;
import com.hartwig.hmftools.svprep.CombinedReadGroups;
import com.hartwig.hmftools.svprep.CombinedStats;
import com.hartwig.hmftools.svprep.ResultsWriter;
import com.hartwig.hmftools.svprep.SvConfig;
import com.hartwig.hmftools.svprep.WriteType;

import htsjdk.samtools.SAMRecord;
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
            final ResultsWriter writer, final CombinedStats combinedStats)
    {
        mId = id;
        mConfig = config;
        mReadFilters = config.ReadFiltering;
        mCombinedReadGroups = combinedReadGroups;
        mWriter = writer;
        mRegion = region;
        mCombinedStats = combinedStats;

        mJunctionTracker = new JunctionTracker(mRegion, mConfig.ReadFiltering.config(), mConfig.Hotspots, mConfig.Blacklist);

        mSamReader = mConfig.BamFile != null ?
                SamReaderFactory.makeDefault().referenceSequence(new File(mConfig.RefGenomeFile)).open(new File(mConfig.BamFile)) : null;

        mBamSlicer = new BamSlicer(0, false, true, false);

        int rateSegmentLength = mConfig.PartitionSize / DOWN_SAMPLE_FRACTION;
        int downsampleThreshold = DOWN_SAMPLE_THRESHOLD / DOWN_SAMPLE_FRACTION;
        mReadRateTracker = new ReadRateTracker(rateSegmentLength, mRegion.start(), downsampleThreshold);
        mRateLimitTriggered = false;

        mFilterRegion = mConfig.RefGenVersion == V37 && region.overlaps(EXCLUDED_REGION_1_REF_37) ? EXCLUDED_REGION_1_REF_37
                : (mConfig.RefGenVersion == V38 && region.overlaps(EXCLUDED_REGION_1_REF_38) ? EXCLUDED_REGION_1_REF_38 : null);

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
                SV_LOGGER.debug("specific readId({})", record.getReadName());
            else if(LOG_READ_ONLY)
                return;
        }

        if(!checkReadRateLimits(readStart))
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
        boolean isSupportCandidate = mReadFilters.isCandidateSupportingRead(record);

        if(!isSupportCandidate && !mConfig.WriteTypes.contains(BAM))
            return;

        ReadRecord read = ReadRecord.from(record);
        read.setFilters(filters);

        if(isSupportCandidate)
            read.setReadType(CANDIDATE_SUPPORT);

        mJunctionTracker.processRead(read);
    }

    private void writeData()
    {
        boolean captureCompleteGroups = mConfig.WriteTypes.contains(BAM) || mConfig.WriteTypes.contains(READS);

        if(captureCompleteGroups)
        {
            // read groups that span chromosomes or partitions need to be complete, so gather up their state to enable this
            Map<String,ReadGroup> spanningGroupsMap = Maps.newHashMap();

            mPerCounters[PC_MATE_SLICE].start();

            mJunctionTracker.junctionGroups().forEach(x -> assignReadGroup(x, spanningGroupsMap));
            mJunctionTracker.supportingGroups().forEach(x -> assignReadGroup(x, spanningGroupsMap));

            mPerCounters[PC_MATE_SLICE].stop();

            int spanningGroups = spanningGroupsMap.size();
            int totalGroups = mJunctionTracker.junctionGroups().size() + mJunctionTracker.supportingGroups().size();

            mCombinedReadGroups.processSpanningReadGroups(mRegion, spanningGroupsMap);

            if(totalGroups == 0)
                return;

            int matchedReads = spanningGroupsMap.values().stream().mapToInt(x -> (int)x.reads().stream().filter(y -> y.written()).count()).sum();

            SV_LOGGER.debug("region({}) readGroups({}) complete(local={} spanning={}) matchedReads({})",
                    mRegion, totalGroups, totalGroups - spanningGroups, spanningGroups, matchedReads);

            if(mConfig.WriteTypes.contains(BAM))
            {
                mWriter.writeBamRecords(mJunctionTracker.junctionGroups());
                mWriter.writeBamRecords(mJunctionTracker.supportingGroups());
            }

            if(mConfig.WriteTypes.contains(READS))
            {
                mWriter.writeReadGroup(mJunctionTracker.junctionGroups());
                mWriter.writeReadGroup(mJunctionTracker.supportingGroups());
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

        if(readGroup.spansPartitions())
        {
            ++mStats.SpanningGroups;

            // get any remote mates
            // findRemoteMateReads(readGroup);
            readGroup.setGroupState();

            // only register remote reads if they are not supplementaries - instead just let these be written locally and not reconciled

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

    private void findRemoteMateReads(final ReadGroup readGroup)
    {
        // find any missing mate reads, excluding supplementaries which are left to be handled individually
        int index = 0;
        while(index < readGroup.reads().size() && index < MAX_GROUP_READ_COUNT)
        {
            ReadRecord read = readGroup.reads().get(index);

            if(!read.isSupplementaryAlignment() && !readGroup.hasReadMate(read))
            {
                SAMRecord mateRecord = mBamSlicer.queryMate(mSamReader, read.record());

                if(mateRecord != null && !readGroup.hasReadRecord(mateRecord))
                {
                    ReadRecord newRead = ReadRecord.from(mateRecord);
                    readGroup.addRead(newRead);
                }
            }

            ++index;
        }
    }

    /*
    private void completeUnmatchedReadGroups(final List<ReadGroupState> unmatchedGroups)
    {
        SV_LOGGER.debug("region({}) extract remote reads for {} unmatched read groups", mRegion, unmatchedGroups.size());

        mPerCounters[PC_UNMATCHED_SLICE].start();

        int completeGroups = 0;

        for(ReadGroupState readGroupState : unmatchedGroups)
        {
            ReadGroup readGroup = mJunctionTracker.getReadGroup(readGroupState.ReadId);

            if(readGroup == null)
            {
                SV_LOGGER.error("readGroup({}) missing, invalid state", readGroupState.ReadId);
                continue;
            }

            completeUnmatchedReadGroup(readGroup);

            if(readGroup.isComplete())
                ++completeGroups;
        }

        SV_LOGGER.debug("region({}) unmatched groups({}) now complete({})", mRegion, unmatchedGroups.size(), completeGroups);

        mPerCounters[PC_UNMATCHED_SLICE].stop();
    }

    private void completeUnmatchedReadGroup(final ReadGroup readGroup)
    {
        // check the need to slice for remote supplementaries
        List<ReadRecord> newReads = Lists.newArrayList();

        for(ReadRecord read : readGroup.reads())
        {
            if(read.hasSuppAlignment())
            {
                SupplementaryReadData suppData = read.supplementaryAlignment();

                if(readGroup.reads().stream().anyMatch(x -> x.Chromosome.equals(suppData.Chromosome) && suppData.Position == read.start()))
                    continue;

                // slice for this supplementary
                GenomePosition slicePosition = GenomePositions.create(suppData.Chromosome, suppData.Position);
                List<SAMRecord> records = mBamSlicer.slice(mSamReader, slicePosition);

                mStats.UnmatchedSliceReads += records.size();

                for(SAMRecord record : records)
                {
                    if(record.getAlignmentStart() != suppData.Position)
                        continue;

                    if(!record.getReadName().equals(readGroup.id()) || !record.getCigarString().equals(suppData.Cigar))
                        continue;

                    ReadRecord newRead = ReadRecord.from(record);
                    int filters = mReadFilters.checkFilters(record);
                    newRead.setFilters(filters);
                    newRead.setReadType(UNMATCHED);
                    newReads.add(newRead);
                    break;
                }
            }
        }

        if(!newReads.isEmpty())
        {
            newReads.forEach(x -> readGroup.addRead(x));
            readGroup.setGroupState();
        }
    }
    */

    private boolean checkReadRateLimits(int positionStart)
    {
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
