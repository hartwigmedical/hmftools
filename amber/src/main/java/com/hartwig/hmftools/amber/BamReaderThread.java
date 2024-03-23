package com.hartwig.hmftools.amber;

import static com.hartwig.hmftools.amber.AmberConfig.AMB_LOGGER;

import java.io.File;
import java.io.IOException;
import java.util.NoSuchElementException;
import java.util.Queue;

import com.hartwig.hmftools.common.bam.BamSlicer;
import com.hartwig.hmftools.common.utils.PerformanceCounter;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;

public class BamReaderThread extends Thread
{
    private final Queue<RegionTask> mTaskQueue;
    private final SamReader mSamReader;
    private final BamSlicer mBamSlicer;
    private RegionTask mCurrentTask;

    private final PerformanceCounter mPerfCounter;

    public BamReaderThread(
            final String bamFile, final SamReaderFactory samReaderFactory, final Queue<RegionTask> inTaskQueue,
            int minMappingQuality)
    {
        mTaskQueue = inTaskQueue;
        mSamReader = samReaderFactory.open(new File(bamFile));
        mBamSlicer = new BamSlicer(minMappingQuality, false, false, false);
        mCurrentTask = null;
        mPerfCounter = new PerformanceCounter("Reads");
    }

    @Override
    public void run()
    {
        // AMB_LOGGER.debug("bam reader thread start");

        while(true)
        {
            RegionTask task;
            try
            {
                task = mTaskQueue.remove();
                mCurrentTask = task;
            }
            catch(NoSuchElementException e)
            {
                // finished processing
                break;
            }

            mPerfCounter.start();
            mBamSlicer.slice(mSamReader, task.Region, this::processRecord);
            mPerfCounter.stop();
        }

        try
        {
            mSamReader.close();
        }
        catch(IOException e)
        {
            AMB_LOGGER.error("IO exception in SamReader::close: {}", e.getMessage());
        }

        // AMB_LOGGER.debug("bam reader thread finish");
    }

    private void processRecord(final SAMRecord record)
    {
        if(mCurrentTask == null)
        {
            mBamSlicer.haltProcessing();
            return;
        }

        mCurrentTask.processRecord(record);

        if(mCurrentTask.isComplete())
            mBamSlicer.haltProcessing();
    }

    public PerformanceCounter perfCounter() { return mPerfCounter; }
}
