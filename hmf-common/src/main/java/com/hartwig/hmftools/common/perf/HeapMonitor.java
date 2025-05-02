package com.hartwig.hmftools.common.perf;

import static java.lang.Math.max;
import static java.lang.String.format;

import java.time.Duration;
import java.util.concurrent.atomic.AtomicBoolean;

import org.apache.logging.log4j.Logger;

public class HeapMonitor extends Thread
{
    private static long MEM_UNIT = 1 << 10;

    private final AtomicBoolean mStopRunning = new AtomicBoolean();

    private final Logger mLogger;
    private final Duration mSamplePeriod;

    private long mMaxUsedMem = 0L;

    public HeapMonitor(final Logger logger, final Duration samplePeriod)
    {
        mLogger = logger;
        mSamplePeriod = samplePeriod;
    }

    private void logMem(final long mem)
    {
        int unitCount = 0;
        float normalisedMem = mem;
        while(normalisedMem >= MEM_UNIT && unitCount < 3)
        {
            unitCount++;
            normalisedMem /= MEM_UNIT;
        }

        String suffix = null;
        if(unitCount == 0)
        {
            suffix = "B";
        }
        else if(unitCount == 1)
        {
            suffix = "KiB";
        }
        else if(unitCount == 2)
        {
            suffix = "MiB";
        }
        else if(unitCount == 3)
        {
            suffix = "GiB";
        }

        mLogger.info(format("Heap usage: %.2f %s", normalisedMem, suffix));
    }

    @Override
    public void run()
    {
        try
        {
            while(!mStopRunning.get())
            {
                System.gc();
                sleep(mSamplePeriod.toMillis());

                Runtime runtime = Runtime.getRuntime();
                long totalMem = runtime.totalMemory();
                long freeMem = runtime.freeMemory();
                long usedMem = totalMem - freeMem;
                mMaxUsedMem = max(mMaxUsedMem, usedMem);
            }

            logMem(mMaxUsedMem);
        }
        catch(InterruptedException e)
        {
            throw new RuntimeException(e);
        }
    }

    public void finish()
    {
        mStopRunning.set(true);
        try
        {
            join();
        }
        catch(InterruptedException e)
        {
            throw new RuntimeException(e);
        }
    }
}
