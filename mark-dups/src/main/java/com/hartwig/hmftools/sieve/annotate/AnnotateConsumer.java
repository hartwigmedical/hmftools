package com.hartwig.hmftools.sieve.annotate;

import static com.hartwig.hmftools.sieve.annotate.AnnotateConfig.MD_LOGGER;

import java.io.File;
import java.util.concurrent.ArrayBlockingQueue;
import java.util.concurrent.Callable;

import com.hartwig.hmftools.common.samtools.BamSlicer;
import com.hartwig.hmftools.common.utils.PerformanceCounter;

import org.jetbrains.annotations.NotNull;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;

public class AnnotateConsumer implements Callable
{
    private final AnnotateConfig mConfig;
    private final ArrayBlockingQueue<IJobRegion> mJobs;

    private final SamReader mSamReader;
    private final BamSlicer mBamSlicer;

    private final PerformanceCounter mConsumerPerfCounter;
    private final PerformanceCounter mJobPerfCounter;

    private IJobRegion mCurrentRegion;
    private long mReadCounter;
    private long mRegionCounter;

    public AnnotateConsumer(@NotNull final AnnotateConfig config, @NotNull final ArrayBlockingQueue<IJobRegion> jobs)
    {
        mConfig = config;
        mJobs = jobs;

        mBamSlicer = new BamSlicer(0, false, true, false);
        mSamReader = SamReaderFactory.makeDefault().referenceSequence(new File(mConfig.RefGenome)).open(new File(mConfig.BamFile));

        mJobPerfCounter = new PerformanceCounter("Job");
        mConsumerPerfCounter = new PerformanceCounter("Total");
        mReadCounter = 0;
        mRegionCounter = 0;
    }

    @Override
    public Long call()
    {
        MD_LOGGER.info("Consumer is starting starting.");

        mConsumerPerfCounter.start();

        while((mCurrentRegion = mJobs.poll()) != null)
        {
            ++mRegionCounter;
            mJobPerfCounter.start();
            mBamSlicer.slice(mSamReader, mCurrentRegion.getChrBaseRegion(), this::processSamRecord);
            mJobPerfCounter.stop();
        }

        mConsumerPerfCounter.stop();

        mJobPerfCounter.logStats();
        mConsumerPerfCounter.logStats();
        MD_LOGGER.info("Consumer is finished, {} reads processed, {} regions processed", mReadCounter, mRegionCounter);

        return (long) 0;
    }

    private void processSamRecord(@NotNull final SAMRecord read)
    {
        ++mReadCounter;

        if(read.getReadUnmappedFlag())
        {
            return;
        }

        mCurrentRegion.matchedRead(read);
    }
}
