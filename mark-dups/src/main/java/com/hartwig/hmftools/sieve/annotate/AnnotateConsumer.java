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
    private final ArrayBlockingQueue<AnnotatedBedRecord> mJobs;

    private final SamReader mSamReader;
    private final BamSlicer mBamSlicer;

    private final PerformanceCounter mPerfCounter;

    private AnnotatedBedRecord mCurrentBedRecord;
    private long mReadCounter;
    private long mBedRecordCounter;

    public AnnotateConsumer(final int taskID, @NotNull final AnnotateConfig config,
            @NotNull final ArrayBlockingQueue<AnnotatedBedRecord> jobs)
    {
        mTaskID = taskID;
        mConfig = config;
        mJobs = jobs;

        mBamSlicer = new BamSlicer(0, false, true, false);
        mSamReader = SamReaderFactory.makeDefault().referenceSequence(new File(mConfig.RefGenome)).open(new File(mConfig.BamFile));

        mPerfCounter = new PerformanceCounter(String.format("Task %d", taskID));
        mReadCounter = 0;
        mBedRecordCounter = 0;
    }

    @Override
    public Long call()
    {
        MD_LOGGER.info("Task {} is starting.", mTaskID);

        AnnotatedBedRecord bedRecord;
        while((bedRecord = mJobs.poll()) != null)
        {
            ++mBedRecordCounter;
            mCurrentBedRecord = bedRecord;

            mPerfCounter.start();

            // TODO(m_cooper): Chunk size.
            ChrBaseRegion partition = new ChrBaseRegion(bedRecord.getChromosome(), bedRecord.getPosStart(), bedRecord.getPosEnd());
            mBamSlicer.slice(mSamReader, partition, this::processSamRecord);

            mPerfCounter.stop();
        }

        mPerfCounter.logStats();
        MD_LOGGER.info("Task {} is finished, {} reads processed, {} BED records processed", mTaskID, mReadCounter, mBedRecordCounter);

        return (long) 0;
    }

    private void processSamRecord(@NotNull final SAMRecord read)
    {
        ++mReadCounter;

        if(read.getReadUnmappedFlag())
        {
            return;
        }

        if(read.getAlignmentStart() >= mCurrentBedRecord.getPosStart() && read.getAlignmentEnd() <= mCurrentBedRecord.getPosEnd())
        {
            mCurrentBedRecord.matchedRead(read);
        }
    }
}
