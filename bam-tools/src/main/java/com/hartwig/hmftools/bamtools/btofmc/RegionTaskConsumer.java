package com.hartwig.hmftools.bamtools.btofmc;

import static com.hartwig.hmftools.bamtools.btofmc.BamToFastqConfig.BFQ_LOGGER;
import static com.hartwig.hmftools.common.bam.SamRecordUtils.CONSENSUS_READ_ATTRIBUTE;

import java.io.IOException;
import java.util.Queue;
import java.util.concurrent.atomic.AtomicInteger;
import java.util.function.Consumer;
import java.util.function.Supplier;

import com.google.common.annotations.VisibleForTesting;
import com.hartwig.hmftools.common.bam.BamSlicer;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReader;

public class RegionTaskConsumer extends Thread
{
    private final BamToFastqConfig mConfig;
    private final Supplier<SamReader> mSamReaderSupplier;
    private final Consumer<SAMRecord> mReadCache;
    private final Queue<RegionTask> mRegionTasks;
    private final AtomicInteger mTotalRegionsProcessedCount;

    private RegionTask mCurrentRegionTask;

    private long mReadsProcessedCount;

    public RegionTaskConsumer(final BamToFastqConfig config, final Supplier<SamReader> samReaderSupplier,
            final Consumer<SAMRecord> readCache, final Queue<RegionTask> regionTasks, final AtomicInteger totalRegionsProcessedCount)
    {
        mConfig = config;
        mSamReaderSupplier = samReaderSupplier;
        mReadCache = readCache;
        mRegionTasks = regionTasks;
        mTotalRegionsProcessedCount = totalRegionsProcessedCount;

        mCurrentRegionTask = null;
        mReadsProcessedCount = 0;
    }

    @VisibleForTesting
    public RegionTaskConsumer(long readsProcessedCount)
    {
        mConfig = null;
        mSamReaderSupplier = null;
        mReadCache = null;
        mRegionTasks = null;
        mTotalRegionsProcessedCount = null;

        mCurrentRegionTask = null;
        mReadsProcessedCount = readsProcessedCount;
    }

    public void run()
    {
        try(SamReader samReader = mSamReaderSupplier.get())
        {
            BamSlicer bamSlicer = new BamSlicer(0, true, false, false);
            bamSlicer.setKeepUnmapped();
            while((mCurrentRegionTask = mRegionTasks.poll()) != null)
            {
                mCurrentRegionTask.slice(samReader, bamSlicer, this::processRead);
                mTotalRegionsProcessedCount.getAndIncrement();
            }
        }
        catch(IOException e)
        {
            BFQ_LOGGER.warn("Failed to close SamReader: {}", e.toString());
        }
    }

    private void processRead(final SAMRecord read)
    {
        if(!mCurrentRegionTask.isAlignmentStartWithinRegion(read))
        {
            return;
        }

        if(!mConfig.KeepConsensusReads && read.hasAttribute(CONSENSUS_READ_ATTRIBUTE))
        {
            return;
        }

        ++mReadsProcessedCount;
        mReadCache.accept(read);
    }

    public void mergeStats(final RegionTaskConsumer other)
    {
        mReadsProcessedCount += other.mReadsProcessedCount;
    }

    public void logStats()
    {
        BFQ_LOGGER.info("RegionTaskConsumer stats: readsProcessedCount({})", mReadsProcessedCount);
    }

    @VisibleForTesting
    public long readsProcessedCount()
    {
        return mReadsProcessedCount;
    }
}
