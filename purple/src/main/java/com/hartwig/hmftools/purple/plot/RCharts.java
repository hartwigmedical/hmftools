package com.hartwig.hmftools.purple.plot;

import java.util.List;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Future;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.utils.r.RExecutor;
import com.hartwig.hmftools.purple.config.ChartConfig;
import com.hartwig.hmftools.purple.config.CommonConfig;
import com.hartwig.hmftools.purple.config.ConfigSupplier;

import org.jetbrains.annotations.NotNull;

class RCharts {

    private final ConfigSupplier configSupplier;
    private final CommonConfig commonConfig;
    private final ChartConfig chartConfig;
    private final ExecutorService executorService;

    RCharts(final ConfigSupplier configSupplier, final ExecutorService executorService) {
        this.configSupplier = configSupplier;
        this.commonConfig = configSupplier.commonConfig();
        this.chartConfig = configSupplier.chartConfig();
        this.executorService = executorService;
    }

    @NotNull
    List<Future<Integer>> chartFutures() {
        final List<Future<Integer>> result = Lists.newArrayList();

        result.add(executorService.submit(() -> RExecutor.executeFromClasspath("r/copyNumberPlots.R",
                commonConfig.tumorSample(),
                commonConfig.outputDirectory(),
                chartConfig.plotDirectory())));

        if (configSupplier.somaticConfig().file().isPresent()) {
            result.add(executorService.submit(() -> RExecutor.executeFromClasspath("r/somaticVariantPlots.R",
                    commonConfig.tumorSample(),
                    commonConfig.outputDirectory(),
                    chartConfig.plotDirectory())));
        }

        return result;
    }
}
