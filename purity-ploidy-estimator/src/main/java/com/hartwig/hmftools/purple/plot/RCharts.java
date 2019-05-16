package com.hartwig.hmftools.purple.plot;

import java.io.IOException;

import com.hartwig.hmftools.common.r.RExecutor;
import com.hartwig.hmftools.purple.config.ChartConfig;
import com.hartwig.hmftools.purple.config.CommonConfig;
import com.hartwig.hmftools.purple.config.ConfigSupplier;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

public class RCharts {

    private static final Logger LOGGER = LogManager.getLogger(RCharts.class);

    private final ConfigSupplier configSupplier;
    private final CommonConfig commonConfig;
    private final ChartConfig chartConfig;

    RCharts(final ConfigSupplier configSupplier) {
        this.configSupplier = configSupplier;
        this.commonConfig = configSupplier.commonConfig();
        this.chartConfig = configSupplier.chartConfig();
    }

    void generateRCharts() throws InterruptedException, IOException {

        int copyNumberResult = RExecutor.executeFromClasspath("r/copyNumberPlots.R",
                commonConfig.tumorSample(),
                commonConfig.outputDirectory(),
                chartConfig.plotDirectory());
        if (copyNumberResult != 0) {
            LOGGER.warn("Error generating R copy number plots.");
        }

        if(configSupplier.somaticConfig().file().isPresent()) {
            int somaticResult = RExecutor.executeFromClasspath("r/somaticVariantPlots.R",
                    commonConfig.tumorSample(),
                    commonConfig.outputDirectory(),
                    chartConfig.plotDirectory());
            if (somaticResult != 0) {
                LOGGER.warn("Error generating R somatic variant plots.");
            }
        }
    }

}
