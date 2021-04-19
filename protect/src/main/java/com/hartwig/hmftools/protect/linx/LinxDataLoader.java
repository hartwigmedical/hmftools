package com.hartwig.hmftools.protect.linx;

import java.io.File;
import java.io.IOException;
import java.util.List;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.variant.structural.linx.LinxBreakend;
import com.hartwig.hmftools.common.variant.structural.linx.LinxFusion;
import com.hartwig.hmftools.common.variant.structural.linx.LinxViralInsertion;
import com.hartwig.hmftools.protect.ProtectConfig;
import com.hartwig.hmftools.protect.purple.PurpleDataLoader;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public final class LinxDataLoader {

    private static final Logger LOGGER = LogManager.getLogger(PurpleDataLoader.class);

    private LinxDataLoader() {
    }

    @NotNull
    public static LinxData load(ProtectConfig config) throws IOException {
        return load(config.linxFusionTsv(), config.linxBreakendTsv(), null, config.linxDriverCatalogTsv());
    }

    @NotNull
    public static LinxData load(@NotNull String linxFusionTsv, @NotNull String linxBreakendTsv, @Nullable String linxViralInsertionTsv,
            @NotNull String linxDriversTsv) throws IOException {
        LOGGER.info("Loading LINX data from {}", new File(linxFusionTsv).getParent());
        List<LinxFusion> linxFusions = LinxFusion.read(linxFusionTsv).stream().filter(LinxFusion::reported).collect(Collectors.toList());
        LOGGER.info(" Loaded {} reportable fusions from {}", linxFusions.size(), linxFusionTsv);

        List<LinxBreakend> linxBreakends =
                LinxBreakend.read(linxBreakendTsv).stream().filter(LinxBreakend::reportedDisruption).collect(Collectors.toList());
        List<ReportableGeneDisruption> reportableGeneDisruptions = ReportableGeneDisruptionFactory.convert(linxBreakends);
        LOGGER.info(" Loaded {} reportable disruptions from {}", reportableGeneDisruptions.size(), linxBreakendTsv);

        // viral insertion is nullable in protect, and notnull in patient reporter
        List<ViralInsertion> reportableViralInsertions = Lists.newArrayList();
        if (linxViralInsertionTsv != null) {
            List<LinxViralInsertion> viralInsertionList = LinxViralInsertion.read(linxViralInsertionTsv);
            reportableViralInsertions = ViralInsertionAnalyzer.analyzeViralInsertions(viralInsertionList);
        }
        LOGGER.info(" Loaded {} reportable viral insertions from {}", reportableViralInsertions.size(), linxViralInsertionTsv);

        List<ReportableHomozygousDisruption> reportableHomozygousDisruptions =
                ReportableHomozygousDisruptionFactory.extractFromLinxDriversTsv(linxDriversTsv);
        LOGGER.info(" Loaded {} reportable homozygous disruptions from {}", reportableHomozygousDisruptions.size(), linxDriversTsv);

        return ImmutableLinxData.builder()
                .addAllFusions(linxFusions)
                .addAllGeneDisruptions(reportableGeneDisruptions)
                .addAllViralInsertions(reportableViralInsertions)
                .addAllHomozygousDisruptions(reportableHomozygousDisruptions)
                .build();
    }
}
