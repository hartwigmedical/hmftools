package com.hartwig.hmftools.sieve.bedstats;

import static java.lang.Math.max;
import static java.lang.Math.min;

import static com.hartwig.hmftools.sieve.annotate.AnnotateConfig.MD_LOGGER;

import java.io.BufferedWriter;
import java.io.File;
import java.util.concurrent.ArrayBlockingQueue;
import java.util.concurrent.Callable;

import com.hartwig.hmftools.common.region.ChrBaseRegion;
import com.hartwig.hmftools.common.samtools.BamSlicer;
import com.hartwig.hmftools.common.utils.PerformanceCounter;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;

public class BedStatsConsumer implements Callable
{
    private final BedStatsConfig mConfig;
    private final ArrayBlockingQueue<ChrBaseRegion> mJobs;
    private final BufferedWriter mOutputWriter;

    private final SamReader mSamReader;
    private final BamSlicer mBamSlicer;

    private final PerformanceCounter mConsumerPerfCounter;
    private final PerformanceCounter mJobPerfCounter;

    private ChrBaseRegion mCurrentRegion;
    private int[] mBaseDepth;
    private long mReadCounter;
    private long mRegionCounter;

    public BedStatsConsumer(final BedStatsConfig config, final ArrayBlockingQueue<ChrBaseRegion> jobs, final BufferedWriter outputWriter)
    {
        mConfig = config;
        mJobs = jobs;
        mOutputWriter = outputWriter;

        mBamSlicer = new BamSlicer(0, false, true, false);
        mSamReader = SamReaderFactory.makeDefault().referenceSequence(new File(mConfig.RefGenome)).open(new File(mConfig.BamFile));

        mBaseDepth = null;
        mJobPerfCounter = new PerformanceCounter("Job");
        mConsumerPerfCounter = new PerformanceCounter("Total");
        mReadCounter = 0;
        mRegionCounter = 0;
    }

    @Override
    public Long call()
    {
        MD_LOGGER.info("Consumer is starting.");

        mConsumerPerfCounter.start();

        while((mCurrentRegion = mJobs.poll()) != null)
        {
            ++mRegionCounter;
            mBaseDepth = new int[mCurrentRegion.baseLength()];

            mJobPerfCounter.start();
            mBamSlicer.slice(mSamReader, mCurrentRegion, this::processSamRecord);
            mJobPerfCounter.stop();

            BedStats.writeRecord(mOutputWriter, mCurrentRegion, new Stats(mBaseDepth));
        }

        mConsumerPerfCounter.stop();

        mJobPerfCounter.logStats();
        mConsumerPerfCounter.logStats();
        MD_LOGGER.info("Consumer is finished, {} reads processed, {} regions processed", mReadCounter, mRegionCounter);

        return (long) 0;
    }

    private void processSamRecord(final SAMRecord read)
    {
        ++mReadCounter;

        if(read.getReadUnmappedFlag())
        {
            return;
        }

        final int readStart = read.getAlignmentStart();
        final int readEnd = read.getAlignmentEnd();
        final int baseStart = max(readStart - mCurrentRegion.start(), 0);
        final int baseEnd = min(readEnd - mCurrentRegion.start(), mBaseDepth.length - 1);

        for(int i = baseStart; i <= baseEnd; ++i)
        {
            ++mBaseDepth[i];
        }
    }
}
