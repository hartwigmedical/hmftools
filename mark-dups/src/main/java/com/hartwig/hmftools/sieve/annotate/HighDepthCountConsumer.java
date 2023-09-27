package com.hartwig.hmftools.sieve.annotate;

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

public class HighDepthCountConsumer implements Callable
{
    private final AnnotateConfig mConfig;
    private final ArrayBlockingQueue<HighDepthRegion> mJobs;
    private final BufferedWriter mOutputWriter;

    private final SamReader mSamReader;
    private final BamSlicer mBamSlicer;

    private final PerformanceCounter mConsumerPerfCounter;
    private final PerformanceCounter mJobPerfCounter;

    private HighDepthRegion mCurrentRegion;
    private HighDepthCounts mCurrentCounts;
    private long mReadCounter;
    private long mRegionCounter;

    public HighDepthCountConsumer(final AnnotateConfig config, final ArrayBlockingQueue<HighDepthRegion> jobs,
            final BufferedWriter outputWriter)
    {
        mConfig = config;
        mJobs = jobs;
        mOutputWriter = outputWriter;

        mBamSlicer = new BamSlicer(0, false, true, false);
        mSamReader = SamReaderFactory.makeDefault().referenceSequence(new File(mConfig.RefGenome)).open(new File(mConfig.BamFile));

        mJobPerfCounter = new PerformanceCounter("HighDepthCountConsumer Jobs");
        mConsumerPerfCounter = new PerformanceCounter("HighDepthCountConsumer Total");
        mReadCounter = 0;
        mRegionCounter = 0;
    }

    @Override
    public Long call()
    {
        mConsumerPerfCounter.start();

        while((mCurrentRegion = mJobs.poll()) != null)
        {
            ++mRegionCounter;
            mCurrentCounts = new HighDepthCounts();
            final ChrBaseRegion chrBaseRegion =
                    new ChrBaseRegion(mCurrentRegion.getChromosome(), mCurrentRegion.getPosStart(), mCurrentRegion.getPosEnd());

            mJobPerfCounter.start();
            mBamSlicer.slice(mSamReader, chrBaseRegion, this::processSamRecord);
            mJobPerfCounter.stop();

            Annotate.writeRecord(mOutputWriter, mCurrentRegion, mCurrentCounts);
        }

        mConsumerPerfCounter.stop();

        mJobPerfCounter.logStats();
        mConsumerPerfCounter.logStats();
        MD_LOGGER.info("HighDepthCountConsumer is finished, {} reads processed, {} regions processed", mReadCounter, mRegionCounter);

        return (long) 0;
    }

    private void processSamRecord(final SAMRecord read)
    {
        ++mReadCounter;

        if(read.getReadUnmappedFlag())
        {
            return;
        }

        mCurrentCounts.matchedRead(read);
    }
}
