package com.hartwig.hmftools.esvee.assembly.types;

import static java.lang.String.format;

import static com.hartwig.hmftools.esvee.AssemblyConfig.SV_LOGGER;

import java.util.List;

import com.hartwig.hmftools.common.utils.PerformanceCounter;

public class ThreadTask extends Thread
{
    protected final PerformanceCounter mPerfCounter;

    public ThreadTask(final String taskName)
    {
        mPerfCounter = new PerformanceCounter(taskName);
    }

    public PerformanceCounter getPerfCounter() { return mPerfCounter; }

    public void stopCheckLog(final String eventDetails, final double threshold)
    {
        mPerfCounter.stop();

        if(threshold > 0 && mPerfCounter.getLastTime() > threshold)
        {
            SV_LOGGER.info("task({}) exceeds perf log time({})", eventDetails, format("%.3f", mPerfCounter.getLastTime()));
        }
    }

    public static void mergePerfCounters(final List<PerformanceCounter> perfCounters, final List<ThreadTask> tasks)
    {
        if(tasks.isEmpty())
            return;

        perfCounters.add(ThreadTask.mergePerfCounters(tasks));
    }

    private static PerformanceCounter mergePerfCounters(final List<ThreadTask> tasks)
    {
        PerformanceCounter first = tasks.get(0).getPerfCounter();

        if(tasks.size() == 1)
            return first;

        PerformanceCounter combined = new PerformanceCounter(first.getName());

        tasks.forEach(x -> combined.merge(x.getPerfCounter()));
        return combined;
    }

}
