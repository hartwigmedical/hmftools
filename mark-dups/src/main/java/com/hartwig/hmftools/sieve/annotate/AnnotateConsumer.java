package com.hartwig.hmftools.sieve.annotate;

import static com.hartwig.hmftools.sieve.annotate.AnnotateConfig.MD_LOGGER;

import java.io.File;
import java.util.concurrent.ArrayBlockingQueue;
import java.util.concurrent.Callable;

import com.hartwig.hmftools.common.region.ChrBaseRegion;
import com.hartwig.hmftools.common.samtools.BamSlicer;
import com.hartwig.hmftools.common.utils.PerformanceCounter;

import org.jetbrains.annotations.NotNull;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;

public class AnnotateConsumer implements Callable
{
    private final int mTaskID;
    private final AnnotateConfig mConfig;
    private final ArrayBlockingQueue<AnnotatedBlacklistRegion> mJobs;

    private final SamReader mSamReader;
    private final BamSlicer mBamSlicer;

    private final PerformanceCounter mPerfCounter;

    private AnnotatedBlacklistRegion mCurrentBlacklistRegion;
    private long mReadCounter;
    private long mBlacklistRegionCounter;

    public AnnotateConsumer(final int taskID, @NotNull final AnnotateConfig config,
            @NotNull final ArrayBlockingQueue<AnnotatedBlacklistRegion> jobs)
    {
        mTaskID = taskID;
        mConfig = config;
        mJobs = jobs;

        mBamSlicer = new BamSlicer(0, false, true, false);
        mSamReader = SamReaderFactory.makeDefault().referenceSequence(new File(mConfig.RefGenome)).open(new File(mConfig.BamFile));

        mPerfCounter = new PerformanceCounter(String.format("Task %d", taskID));
        mReadCounter = 0;
        mBlacklistRegionCounter = 0;
    }

    @Override
    public Long call()
    {
        MD_LOGGER.info("Consumer is starting starting.");

        AnnotatedBlacklistRegion blacklistRegion;
        while((blacklistRegion = mJobs.poll()) != null)
        {
            ++mBlacklistRegionCounter;
            mCurrentBlacklistRegion = blacklistRegion;

            mPerfCounter.start();

            // TODO(m_cooper): Chunk size.
            ChrBaseRegion partition = new ChrBaseRegion(
                    mCurrentBlacklistRegion.getBlacklistRegion().getChromosome(),
                    mCurrentBlacklistRegion.getBlacklistRegion().getPosStart(),
                    mCurrentBlacklistRegion.getBlacklistRegion().getPosEnd());
            mBamSlicer.slice(mSamReader, partition, this::processSamRecord);

            mPerfCounter.stop();
        }

        mPerfCounter.logStats();
        MD_LOGGER.info("Consumer is finished, {} reads processed, {} blacklisted regions processed", mReadCounter, mBlacklistRegionCounter);

        return (long) 0;
    }

    private void processSamRecord(@NotNull final SAMRecord read)
    {
        ++mReadCounter;

        if(read.getReadUnmappedFlag())
        {
            return;
        }

        if(read.getAlignmentStart() >= mCurrentBlacklistRegion.getBlacklistRegion().getPosStart()
                && read.getAlignmentEnd() <= mCurrentBlacklistRegion.getBlacklistRegion().getPosEnd())
        {
            mCurrentBlacklistRegion.matchedRead(read);
        }
    }
}
