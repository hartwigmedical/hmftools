package com.hartwig.hmftools.svtools.sv_prep;

import static com.hartwig.hmftools.svtools.sv_prep.SvCommon.SV_LOGGER;

import java.io.File;
import java.util.Map;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.samtools.BamSlicer;
import com.hartwig.hmftools.common.samtools.SupplementaryReadData;
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

    private final SvBucket[] mBuckets;
    private int mProcessedBucketIndex;

    private final PartitionStats mStats;

    private final Map<String,ReadGroup> mReadGroups;

    public PartitionSlicer(final int id, final ChrBaseRegion region, final SvConfig config, final ResultsWriter writer)
    {
        mId = id;
        mConfig = config;
        mReadFilters = config.ReadFilters;
        mWriter = writer;
        mRegion = region;

        mSamReader = mConfig.BamFile != null ?
                SamReaderFactory.makeDefault().referenceSequence(new File(mConfig.RefGenomeFile)).open(new File(mConfig.BamFile)) : null;

        mBamSlicer = new BamSlicer(0, false, true, false);

        int bucketCount = mConfig.PartitionSize / mConfig.BucketSize + 1;
        mBuckets = new SvBucket[bucketCount];
        mProcessedBucketIndex = -1;

        mReadGroups = Maps.newHashMap();

        mStats = new PartitionStats();
    }

    public void run()
    {
        SV_LOGGER.debug("processing region({})", mRegion);

        mBamSlicer.slice(mSamReader, Lists.newArrayList(mRegion), this::processSamRecord);

        mReadGroups.values().forEach(x -> processGroup(x));

        processBuckets(-1);

        SV_LOGGER.debug("region({}) complete, stats({}) filters({}) incompleteGroups({})",
                mRegion, mStats.toString(), ReadFilterType.filterCountsToString(mStats.ReadFilterCounts), mReadGroups.size());
    }

    private void processSamRecord(final SAMRecord record)
    {
        if(!mRegion.containsPosition(record.getAlignmentStart()))
            return;

        ++mStats.TotalReads;

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
            mReadGroups.put(read.Id, new ReadGroup(read));
            return;
        }

        readGroup.addRead(read);

        boolean processGroup = readGroup.isComplete();

        if(!processGroup)
        {
            // check if the supplement is also in this partition
            SupplementaryReadData suppData = read.supplementaryAlignment();

            if(suppData == null || !mRegion.containsPosition(suppData.Chromosome, suppData.Position))
                processGroup = true;
        }

        if(processGroup)
            processGroup(readGroup);

        if(readGroup.isComplete())
            mReadGroups.remove(readGroup.id());
    }

    private void processSingleRead(final ReadRecord read)
    {
        SvBucket bucket = findBucket(read.start());
        bucket.addReadGroup(new ReadGroup(read));
    }

    private void processGroup(final ReadGroup readGroup)
    {
        SvBucket bucket = findBucket(readGroup.minStartPosition());
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
        SvBucket bucket = findBucket(read.start());
        bucket.supportingReads().add(read);
    }

    private void processBuckets(int currentReadPosition)
    {
        while(mProcessedBucketIndex < mBuckets.length - 1)
        {
            int nextBucketIndex = mProcessedBucketIndex + 1;
            int bucketStartPos = nextBucketIndex * mConfig.BucketSize;

            if(currentReadPosition > 0 && currentReadPosition < bucketStartPos - mConfig.BucketSize)
                break;

            mProcessedBucketIndex = nextBucketIndex;
            SvBucket bucket = mBuckets[mProcessedBucketIndex];

            if(bucket == null)
                continue;

            processBucket(bucket);
        }
    }

    private void processBucket(final SvBucket bucket)
    {
        if(bucket.junctionPositions().isEmpty())
        {
            ++mStats.EmptyBuckets;
            return;
        }

        ++mStats.Buckets;

        bucket.setJunctionPositions();
        bucket.selectSupportingReads();

        mStats.JunctionCount += bucket.junctionPositions().size();
        mStats.JunctionFragmentCount += bucket.readGroups().size();
        mStats.SupportingReadCount += bucket.supportingReads().size();
        mStats.InitialSupportingReadCount = bucket.initialSupportingReadCount();

        mWriter.writeBucketData(bucket, mId);

        if(mConfig.WriteTypes.contains(WriteType.READS))
        {
            bucket.readGroups().forEach(x -> x.reads().forEach(y -> mWriter.writeReadData(y, "SV", mId, bucket.id())));
            bucket.supportingReads().forEach(x -> mWriter.writeReadData(x, "SUPPORTING", mId, bucket.id()));
        }
    }

    private SvBucket findBucket(int position)
    {
        int positionOffset = position - mRegion.start();

        int bucket = positionOffset / mConfig.BucketSize;

        if(bucket < 0 || bucket >= mBuckets.length)
        {
            SV_LOGGER.error("partition({}) invalid bucket index({}) from position({} offset={}) vs max({})",
                    mRegion, bucket, position, positionOffset, mBuckets.length);
            return null;
        }

        if(mBuckets[bucket] == null)
        {
            mBuckets[bucket] = new SvBucket(bucket);
        }

        return mBuckets[bucket];
    }
}
