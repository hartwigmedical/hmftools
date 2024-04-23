package com.hartwig.hmftools.bamtools.tofastq;

import static java.lang.String.format;

import static com.hartwig.hmftools.bamtools.common.CommonUtils.BT_LOGGER;
import static com.hartwig.hmftools.bamtools.tofastq.PartitionData.formChromosomePartition;
import static com.hartwig.hmftools.common.bam.SamRecordUtils.CONSENSUS_READ_ATTRIBUTE;

import java.io.File;
import java.io.IOException;
import java.util.List;
import java.util.Map;
import java.util.NoSuchElementException;
import java.util.Queue;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.bam.BamSlicer;
import com.hartwig.hmftools.common.region.ChrBaseRegion;
import com.hartwig.hmftools.common.utils.PerformanceCounter;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;

public class PartitionReader extends Thread
{
    private final FastqConfig mConfig;

    private final PartitionDataStore mPartitionDataStore;

    private final FastqWriterCache mWriterCache;
    private final FastqWriter mThreadedWriter;

    private final Queue<ChrBaseRegion> mPartitions;
    private final int mTotalPartitionCount;

    private final SamReader mSamReader;
    private final BamSlicer mBamSlicer;

    private ChrBaseRegion mCurrentRegion;

    private String mCurrentStrPartition;
    private PartitionData mCurrentPartitionData;

    private final Map<String,List<SAMRecord>> mRemoteUnmatchedReads;
    private final Map<String,SAMRecord> mLocalUnmatchedReads;

    private int mPartitionRecordCount;
    private final PerformanceCounter mPerfCounter;

    public PartitionReader(
            final FastqConfig config, final Queue<ChrBaseRegion> partitions,
            final FastqWriterCache writerCache, final PartitionDataStore partitionDataStore)
    {
        mConfig = config;
        mPartitionDataStore = partitionDataStore;

        mWriterCache = writerCache;

        mThreadedWriter = mConfig.SplitMode == FileSplitMode.THREAD ? writerCache.createThreadedWriter() : null;

        mPartitions = partitions;
        mTotalPartitionCount = partitions.size();

        mSamReader = mConfig.BamFile != null ?
                SamReaderFactory.makeDefault().referenceSequence(new File(mConfig.RefGenomeFile)).open(new File(mConfig.BamFile)) : null;

        mBamSlicer = new BamSlicer(0, true, false, false);
        mBamSlicer.setKeepUnmapped();

        mCurrentRegion = null;

        mLocalUnmatchedReads = Maps.newHashMap();
        mRemoteUnmatchedReads = Maps.newHashMap();

        mPartitionRecordCount = 0;

        mPerfCounter = new PerformanceCounter("PartitionReader");
    }

    private static final int LOG_COUNT = 100;

    @Override
    public void run()
    {
        if(mPartitions.isEmpty())
            return;

        while(true)
        {
            try
            {
                ChrBaseRegion region = mPartitions.remove();

                int remainingCount = mPartitions.size();
                int processedCount = mTotalPartitionCount - remainingCount;

                if(processedCount > 0 && (processedCount % LOG_COUNT) == 0)
                {
                    BT_LOGGER.debug("processing partition({}), remaining({})", region, remainingCount);
                }

                processRegion(region);
            }
            catch(NoSuchElementException e)
            {
                BT_LOGGER.trace("all tasks complete");
                break;
            }
            catch(Exception e)
            {
                e.printStackTrace();
                System.exit(1);
            }
        }

        try
        {
            mSamReader.close();
        }
        catch(IOException e)
        {
            BT_LOGGER.error("failed to close bam file: {}", e.toString());
        }
    }

    public PerformanceCounter perfCounter() { return mPerfCounter; }

    public void processRegion(final ChrBaseRegion region)
    {
        mCurrentRegion = region;

        perfCountersStart();

        mCurrentStrPartition = formChromosomePartition(region.Chromosome, mCurrentRegion.start(), mConfig.PartitionSize);
        mCurrentPartitionData = mPartitionDataStore.getOrCreatePartitionData(mCurrentStrPartition);

        if(mBamSlicer != null)
        {
            mBamSlicer.slice(mSamReader, mCurrentRegion, this::processSamRecord);
        }

        postProcessRegion();
    }

    @VisibleForTesting
    public void postProcessRegion()
    {
        // post-slice clean-up
        processPendingIncompletes();

        perfCountersStop();

        BT_LOGGER.debug("partition({}) complete, reads({})", mCurrentRegion, mPartitionRecordCount);

        if(mConfig.PerfDebug)
            mCurrentPartitionData.logCacheCounts();

        mPartitionRecordCount = 0;
    }

    private void processSamRecord(final SAMRecord read)
    {
        if(!mCurrentRegion.containsPosition(read.getAlignmentStart())) // to avoid processing reads from the prior region again
            return;

        if(read.hasAttribute(CONSENSUS_READ_ATTRIBUTE)) // drop any consensus reads
            return;

        if(!read.getReadPairedFlag())
        {
            mWriterCache.processUnpairedRead(read, mThreadedWriter);
            return;
        }

        SAMRecord mate = mLocalUnmatchedReads.remove(read.getReadName());

        if(mate != null)
        {
            mWriterCache.processReadPair(read, mate, mThreadedWriter);
            return;
        }

        boolean hasLocalMate = read.getReferenceIndex() == read.getMateReferenceIndex()
                && mCurrentRegion.containsPosition(read.getMateAlignmentStart());

        if(hasLocalMate)
        {
            mLocalUnmatchedReads.put(read.getReadName(), read);
            return;
        }

        // cache if local otherwise put into remote pending list
        processRemoteRead(read);

        ++mPartitionRecordCount;
    }

    public String getBasePartition(final SAMRecord read)
    {
        // take the lower of the read and its mate
        boolean readLowerPos;
        if(read.getReferenceIndex() == read.getMateReferenceIndex())
        {
            readLowerPos = read.getAlignmentStart() < read.getMateAlignmentStart();
        }
        else
        {
            readLowerPos = read.getReferenceIndex() < read.getMateReferenceIndex();
        }

        return readLowerPos ?
                formChromosomePartition(read.getReferenceName(), read.getAlignmentStart(), mConfig.PartitionSize)
                : formChromosomePartition(read.getMateReferenceName(), read.getMateAlignmentStart(), mConfig.PartitionSize);
    }

    private void processRemoteRead(final SAMRecord read)
    {
        String basePartition = getBasePartition(read);

        // ++mStats.InterPartition;

        // cache this read and send through as groups when the partition is complete
        List<SAMRecord> pendingFragments = mRemoteUnmatchedReads.get(basePartition);

        if(pendingFragments == null)
        {
            pendingFragments = Lists.newArrayList();
            mRemoteUnmatchedReads.put(basePartition, pendingFragments);
        }

        pendingFragments.add(read);
    }

    private void processPendingIncompletes()
    {
        if(mRemoteUnmatchedReads.isEmpty())
            return;

        if(mRemoteUnmatchedReads.size() > 10000)
        {
            BT_LOGGER.debug("partition({}) processing {} remote reads from {} remote partitions",
                    mCurrentRegion, mRemoteUnmatchedReads.values().stream().mapToInt(x -> x.size()).sum(), mRemoteUnmatchedReads.size());
        }

        for(Map.Entry<String,List<SAMRecord>> entry : mRemoteUnmatchedReads.entrySet())
        {
            String basePartition = entry.getKey();
            List<SAMRecord> reads = entry.getValue();

            PartitionData partitionData = mPartitionDataStore.getOrCreatePartitionData(basePartition);

            List<ReadPair> matchedPairs = partitionData.processUnpairedReads(reads);

            for(ReadPair readPair : matchedPairs)
            {
                mWriterCache.processReadPair(readPair.First, readPair.Second, mThreadedWriter);
            }
        }

        mRemoteUnmatchedReads.clear();
    }

    private void perfCountersStart()
    {
        if(mConfig.PerfDebug)
            mPerfCounter.start(format("%s", mCurrentRegion));
        else
            mPerfCounter.start();
    }

    private void perfCountersStop()
    {
        mPerfCounter.stop();
    }

    @VisibleForTesting
    public void processRead(final SAMRecord read) { processSamRecord(read); }

    @VisibleForTesting
    public void flushPendingIncompletes() { processPendingIncompletes(); }

    @VisibleForTesting
    public PartitionDataStore partitionDataStore() { return mPartitionDataStore; }
}
