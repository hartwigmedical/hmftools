package com.hartwig.hmftools.purple.plot;

import java.io.File;
import java.io.IOException;
import java.util.List;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Future;

import com.hartwig.hmftools.common.amber.AmberBAF;
import com.hartwig.hmftools.common.purple.copynumber.PurpleCopyNumber;
import com.hartwig.hmftools.common.purple.gender.Gender;
import com.hartwig.hmftools.common.purple.region.FittedRegion;
import com.hartwig.hmftools.common.variant.structural.StructuralVariant;
import com.hartwig.hmftools.purple.config.ChartConfig;
import com.hartwig.hmftools.purple.config.ConfigSupplier;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

import htsjdk.variant.variantcontext.VariantContext;

public class Charts {

    private static final Logger LOGGER = LogManager.getLogger(Charts.class);

    private final ConfigSupplier configSupplier;
    private final ExecutorService executorService;

    public Charts(final ConfigSupplier configSupplier, final ExecutorService executorService) throws IOException {
        this.configSupplier = configSupplier;
        this.executorService = executorService;

        ChartConfig chartConfig = configSupplier.chartConfig();
        createDirectory(chartConfig.circosDirectory());
        if (chartConfig.enabled() || chartConfig.circosBinary().isPresent()) {
            createDirectory(configSupplier.chartConfig().plotDirectory());
        }
    }

    public void write(@NotNull final Gender gender, @NotNull final List<PurpleCopyNumber> copyNumbers,
            @NotNull final List<VariantContext> somaticVariants, @NotNull final List<StructuralVariant> structuralVariants,
            @NotNull final List<FittedRegion> regions, @NotNull final List<AmberBAF> bafs)
            throws InterruptedException, ExecutionException, IOException {

        final ChartConfig chartConfig = configSupplier.chartConfig();
        final CircosCharts circosCharts = new CircosCharts(configSupplier, executorService);
        circosCharts.write(gender, copyNumbers, somaticVariants, structuralVariants, regions, bafs);
        final List<Future<Integer>> futures = circosCharts.chartFutures();

        if (chartConfig.enabled()) {
            final RCharts rCharts = new RCharts(configSupplier, executorService);
            futures.addAll(rCharts.chartFutures());
        }

        for (final Future<Integer> future : futures) {
            // This (intentionally) has side effect of alerting users to any exceptions
            int result = future.get();
            if (result != 0) {
                LOGGER.warn("Error generating charts.");
            }
        }
    }

    private void createDirectory(final String dir) throws IOException {
        final File output = new File(dir);
        if (!output.exists() && !output.mkdirs()) {
            throw new IOException("Unable to create chart directory " + dir);
        }
    }
}
