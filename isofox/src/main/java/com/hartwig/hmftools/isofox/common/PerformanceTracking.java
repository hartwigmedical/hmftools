package com.hartwig.hmftools.isofox.common;

import static java.lang.Math.round;

import static com.hartwig.hmftools.isofox.IsofoxConfig.ISF_LOGGER;

import java.util.Arrays;
import java.util.List;

import com.hartwig.hmftools.common.utils.PerformanceCounter;
import com.hartwig.hmftools.isofox.IsofoxConfig;

public class PerformanceTracking
{
    private final IsofoxConfig mConfig;

    // suite of performance trackers
    public static final int PERF_TOTAL = 0;
    public static final int PERF_READS = 1;
    public static final int PERF_NOVEL_LOCATIONS = 2;
    public static final int PERF_FIT = 3;
    public static final int PERF_GC_ADJUST = 4;
    public static final int PERF_FUSIONS = 5;

    private static final int PERF_MAX = PERF_FUSIONS+1;

    public PerformanceTracking(final IsofoxConfig config)
    {
        mConfig = config;
    }

    public boolean enabled() { return mConfig.RunPerfChecks; }


    public final static PerformanceCounter[] createPerfCounters()
    {
        PerformanceCounter[] perfCounters = new PerformanceCounter[PERF_MAX];
        perfCounters[PERF_TOTAL] = new PerformanceCounter("Total");
        perfCounters[PERF_READS] = new PerformanceCounter("ReadCounts");
        perfCounters[PERF_NOVEL_LOCATIONS] = new PerformanceCounter("NovelLocations");
        perfCounters[PERF_FIT] = new PerformanceCounter("ExpressFit");
        perfCounters[PERF_FUSIONS] = new PerformanceCounter("Fusions");
        perfCounters[PERF_GC_ADJUST] = new PerformanceCounter("GcAdjust");

        //if(mConfig.RunPerfChecks)
        //    perfCounters[PERF_FIT].setSortTimes(true);

        return perfCounters;
    }

    public void logPerformanceStats(final List<PerformanceCounter[]> perfCounters)
    {
        if(!mConfig.RunPerfChecks)
            return;

        final PerformanceCounter[] combinedPc = perfCounters.get(0);

        for(int i = 1; i < perfCounters.size(); ++i)
        {
            final PerformanceCounter[] chrPCs = perfCounters.get(i);

            for(int j = 0; j < combinedPc.length; ++j)
            {
                combinedPc[j].merge(chrPCs[j]);
            }
        }

        Arrays.stream(combinedPc).forEach(x -> x.logStats());

        /*
        if(mConfig.RunPerfChecks)
        {
            // log 10 slowest times and their interval names
            final List<Double> fitTimes = combinedPc[PERF_FIT].getTimes();
            final List<String> fitGenes = combinedPc[PERF_FIT].getTimeNames();

            if(fitTimes.size() >= 10 && fitGenes.size() == fitTimes.size())
            {
                for (int i = fitTimes.size() - 1; i >= fitTimes.size() - 10; --i)
                {
                    ISF_LOGGER.info(String.format("fit times: geneSet(%s) time(%.3f)", fitGenes.get(i), fitTimes.get(i)));
                }
            }
        }
        */
    }

    private static final long MEGABYTE = 1024L * 1024L;

    public static int calcCurrentMemoryUsage(boolean runGc)
    {
        Runtime runtime = Runtime.getRuntime();

        if(runGc)
            runtime.gc();

        long memory = runtime.totalMemory() - runtime.freeMemory();
        return round(memory / MEGABYTE);
    }

    public static void logMemory(final IsofoxConfig config, final String stage)
    {
        if(!config.RunPerfChecks)
            return;

        ISF_LOGGER.info("stage({}) memory({})", stage, calcCurrentMemoryUsage(true));
    }

}
