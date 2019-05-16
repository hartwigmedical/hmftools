package com.hartwig.hmftools.purple.plot;

import java.io.File;
import java.io.IOException;
import java.util.List;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;

import com.hartwig.hmftools.common.amber.AmberBAF;
import com.hartwig.hmftools.common.purple.copynumber.PurpleCopyNumber;
import com.hartwig.hmftools.common.purple.gender.Gender;
import com.hartwig.hmftools.common.purple.region.FittedRegion;
import com.hartwig.hmftools.common.variant.PurityAdjustedSomaticVariant;
import com.hartwig.hmftools.common.variant.structural.StructuralVariant;
import com.hartwig.hmftools.purple.config.ChartConfig;
import com.hartwig.hmftools.purple.config.ConfigSupplier;

import org.jetbrains.annotations.NotNull;

public class Charts {

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
            @NotNull final List<PurityAdjustedSomaticVariant> somaticVariants, @NotNull final List<StructuralVariant> structuralVariants,
            @NotNull final List<FittedRegion> regions, @NotNull final List<AmberBAF> bafs)
            throws InterruptedException, ExecutionException, IOException {

        final ChartConfig chartConfig = configSupplier.chartConfig();
        new CircosCharts(configSupplier, executorService).write(gender, copyNumbers, somaticVariants, structuralVariants, regions, bafs);

        if (chartConfig.enabled()) {
            new RCharts(configSupplier).generateRCharts();
        }

    }

    private void createDirectory(final String dir) throws IOException {
        final File output = new File(dir);
        if (!output.exists() && !output.mkdirs()) {
            throw new IOException("Unable to create chart directory " + dir);
        }
    }

}
