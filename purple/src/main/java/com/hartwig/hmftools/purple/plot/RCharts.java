package com.hartwig.hmftools.purple.plot;

import java.util.List;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Future;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.utils.r.RExecutor;
import com.hartwig.hmftools.purple.ChartConfig;
import com.hartwig.hmftools.purple.PurpleConfig;

public class RCharts
{
    private final ChartConfig mChartConfig;
    private final ExecutorService mExecutorService;
    private final String mOutputDir;

    public RCharts(final PurpleConfig config, final ExecutorService executorService)
    {
        mChartConfig = config.Charting;
        mExecutorService = executorService;
        mOutputDir = config.OutputDir;
    }

    public List<Future<Integer>> chartFutures(final String sampleId, boolean plotSomatics)
    {
        final List<Future<Integer>> result = Lists.newArrayList();

        result.add(mExecutorService.submit(() -> RExecutor.executeFromClasspath("r/copyNumberPlots.R",
                sampleId,
                mOutputDir,
                mChartConfig.PlotDirectory)));

        if(plotSomatics)
        {
            result.add(mExecutorService.submit(() -> RExecutor.executeFromClasspath("r/somaticVariantPlots.R",
                    sampleId,
                    mOutputDir,
                    mChartConfig.PlotDirectory)));
        }

        return result;
    }
}
