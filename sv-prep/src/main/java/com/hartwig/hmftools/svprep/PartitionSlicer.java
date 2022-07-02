package com.hartwig.hmftools.svprep;

import static com.hartwig.hmftools.svprep.SvCommon.SV_LOGGER;
import static com.hartwig.hmftools.svprep.SvConstants.DOWN_SAMPLE_FRACTION;
import static com.hartwig.hmftools.svprep.SvConstants.DOWN_SAMPLE_THRESHOLD;

import java.io.File;
import java.util.Map;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.samtools.BamSlicer;
import com.hartwig.hmftools.common.samtools.SupplementaryReadData;
import com.hartwig.hmftools.common.utils.PerformanceCounter;
import com.hartwig.hmftools.common.utils.sv.ChrBaseRegion;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;

public class PartitionSlicer
{
    private final int mId;
    private final SvConfig mConfig;
    private final ChrBaseRegion mRegion;
    private final ResultsWriter mWriter;
    private final ReadFilterConfig mReadFilters;

    private final SamReader mSamReader;
    private final BamSlicer mBamSlicer;

    private final PartitionBuckets mBuckets;

    private final PartitionStats mStats;
    private final CombinedStats mCombinedStats;

    private final Map<String,ReadGroup> mReadGroups;

    private final ReadRateTracker mReadRateTracker;
    private boolean mRateLimitTriggered;
    private boolean mLogReadIds;
    private final PerformanceCounter mPerCounter;

    public PartitionSlicer(
            final int id, final ChrBaseRegion region, final SvConfig config, final ResultsWriter writer, final CombinedStats combinedStats)
    {
        mId = id;
        mConfig = config;
        mReadFilters = config.ReadFilters;
        mWriter = writer;
        mRegion = region;
        mCombinedStats = combinedStats;

        mSamReader = mConfig.BamFile != null ?
                SamReaderFactory.makeDefault().referenceSequence(new File(mConfig.RefGenomeFile)).open(new File(mConfig.BamFile)) : null;

        mBamSlicer = new BamSlicer(0, false, true, false);

        mBuckets = new PartitionBuckets(mRegion, mConfig.PartitionSize, mConfig.BucketSize);

        int rateSegmentLength = mConfig.PartitionSize / DOWN_SAMPLE_FRACTION;
        int downsampleThreshold = DOWN_SAMPLE_THRESHOLD / DOWN_SAMPLE_FRACTION;
        mReadRateTracker = new ReadRateTracker(rateSegmentLength, mRegion.start(), downsampleThreshold);
        mRateLimitTriggered = false;

        mReadGroups = Maps.newHashMap();

        mStats = new PartitionStats();
        mPerCounter = new PerformanceCounter("Total");
        mLogReadIds = !mConfig.LogReadIds.isEmpty();
    }

    public void run()
    {
        SV_LOGGER.debug("processing region({})", mRegion);

        mPerCounter.start();

        mBamSlicer.slice(mSamReader, Lists.newArrayList(mRegion), this::processSamRecord);

        mReadGroups.values().forEach(x -> processGroup(x));

        mBuckets.processBuckets(-1, this::processBucket);

        mPerCounter.stop();

        if(mStats.TotalReads == 0)
            return;

        SV_LOGGER.debug("region({}) complete, stats({}) incompleteGroups({})",
                mRegion, mStats.toString(), mReadGroups.size());

        SV_LOGGER.debug("region({}) filters({})",
                mRegion, ReadFilterType.filterCountsToString(mStats.ReadFilterCounts));

        mCombinedStats.addPartitionStats(mStats);
        mCombinedStats.addPerfCounters(mPerCounter);

        if(mRateLimitTriggered)
            System.gc();
    }

    private static final boolean LOG_READ_ONLY = false;

    private void processSamRecord(final SAMRecord record)
    {
        if(!mRegion.containsPosition(record.getAlignmentStart()))
            return;

        ++mStats.TotalReads;

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

        if(!checkReadRateLimits(record.getAlignmentStart()))
            return;

        int filters = mReadFilters.checkFilters(record);

        if(filters != 0)
        {
            processFilteredRead(record, filters);
            return;
        }

        ReadRecord read = ReadRecord.from(record);

        if(!mRegion.containsPosition(read.MateChromosome, read.MatePosStart))
        {
            processSingleRead(read);
            return;
        }

        ReadGroup readGroup = mReadGroups.get(read.Id);

        if(readGroup == null)
        {
            // cache the read waiting for its mate
            mReadGroups.put(read.Id, new ReadGroup(read));
            return;
        }

        readGroup.addRead(read);

        if(readGroup.isComplete())
        {
            processGroup(readGroup);
            mReadGroups.remove(readGroup.id());
            return;
        }

        // if either read has a supplementary in another partition then process this incomplete group now
        SupplementaryReadData suppData = readGroup.reads().stream()
                .filter(x -> x.hasSuppAlignment()).findFirst().map(x -> x.supplementaryAlignment()).orElse(null);

        if(suppData != null && !mRegion.containsPosition(suppData.Chromosome, suppData.Position))
        {
            processGroup(readGroup);
            mReadGroups.remove(readGroup.id());
        }
    }

    private void processSingleRead(final ReadRecord read)
    {
        SvBucket bucket = mBuckets.findBucket(read.start());
        bucket.addReadGroup(new ReadGroup(read));
    }

    private void processGroup(final ReadGroup readGroup)
    {
        SvBucket bucket = mBuckets.findBucket(readGroup.minStartPosition());
        bucket.addReadGroup(readGroup);
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
        if(!mReadFilters.isCandidateSupportingRead(record))
            return;

        ReadRecord read = ReadRecord.from(record);
        SvBucket bucket = mBuckets.findBucket(read.start());
        bucket.addSupportingRead(read);
    }

    private void processBucket(final SvBucket bucket)
    {
        bucket.createJunctions();

        // apply basic filters
        bucket.filterJunctions(mConfig.Hotspots, mReadFilters.MinJunctionSupport);

        if(bucket.junctions().isEmpty())
        {
            ++mStats.FilteredBuckets;
            return;
        }

        ++mStats.Buckets;

        mStats.JunctionCount += bucket.junctions().size();
        mStats.JunctionFragmentCount += bucket.readGroupCount();
        mStats.SupportingReadCount += bucket.supportingReads().size();
        mStats.InitialSupportingReadCount += bucket.initialSupportingReadCount();

        if(mConfig.WriteTypes.contains(WriteType.BUCKET_STATS))
        {
            mWriter.writeBucketData(bucket, mId);
        }

        if(mConfig.WriteTypes.contains(WriteType.JUNCTIONS))
        {
            mWriter.writeJunctionData(bucket, mId);
        }

        if(mConfig.WriteTypes.contains(WriteType.READS))
        {
            for(ReadGroup readGroup : bucket.readGroups())
            {
                boolean groupComplete = readGroup.isComplete();
                readGroup.reads().forEach(x -> mWriter.writeReadData(x, mId, bucket.id(), "SV", groupComplete));
            }

            // group complete set false since reads are not grouped for now
            bucket.supportingReads().forEach(x -> mWriter.writeReadData(x, mId, bucket.id(), "SUPPORTING", false));
        }
    }

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
