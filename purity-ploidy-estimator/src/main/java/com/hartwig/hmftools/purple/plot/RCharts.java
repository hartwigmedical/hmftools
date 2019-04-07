package com.hartwig.hmftools.purple.plot;

import java.io.IOException;

import com.hartwig.hmftools.common.r.RExecutor;
import com.hartwig.hmftools.purple.config.ChartConfig;
import com.hartwig.hmftools.purple.config.CommonConfig;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

public class RCharts {

    private static final Logger LOGGER = LogManager.getLogger(RCharts.class);

    private final CommonConfig commonConfig;
    private final ChartConfig chartConfig;

    public RCharts(final CommonConfig commonConfig, final ChartConfig chartConfig) {
        this.commonConfig = commonConfig;
        this.chartConfig = chartConfig;
    }

    public void generatePlots() throws InterruptedException, IOException {

        int result = RExecutor.executeFromClasspath("r/purplePlots.R",
                commonConfig.tumorSample(),
                commonConfig.outputDirectory(),
                chartConfig.plotDirectory());
        if (result != 0) {
            LOGGER.warn("Error generating R plots.");
        }

    }

}
