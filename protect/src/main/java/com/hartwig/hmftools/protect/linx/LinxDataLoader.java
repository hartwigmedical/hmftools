package com.hartwig.hmftools.protect.linx;

import java.io.File;
import java.io.IOException;
import java.util.List;
import java.util.stream.Collectors;

import com.hartwig.hmftools.common.variant.structural.linx.LinxBreakend;
import com.hartwig.hmftools.common.variant.structural.linx.LinxFusion;
import com.hartwig.hmftools.protect.ProtectConfig;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

public final class LinxDataLoader {

    private static final Logger LOGGER = LogManager.getLogger(LinxDataLoader.class);

    private LinxDataLoader() {
    }

    @NotNull
    public static LinxData load(@NotNull ProtectConfig config) throws IOException {
        return load(config.linxFusionTsv(), config.linxBreakendTsv(), config.linxDriverCatalogTsv());
    }

    @NotNull
    public static LinxData load(@NotNull String linxFusionTsv, @NotNull String linxBreakendTsv,
            @NotNull String linxDriverCatalogTsv) throws IOException {
        LOGGER.info("Loading LINX data from {}", new File(linxFusionTsv).getParent());
        List<LinxFusion> linxFusions = LinxFusion.read(linxFusionTsv).stream().filter(LinxFusion::reported).collect(Collectors.toList());
        LOGGER.info(" Loaded {} reportable fusions from {}", linxFusions.size(), linxFusionTsv);

        List<LinxBreakend> linxBreakends =
                LinxBreakend.read(linxBreakendTsv).stream().filter(LinxBreakend::reportedDisruption).collect(Collectors.toList());
        List<ReportableGeneDisruption> reportableGeneDisruptions = ReportableGeneDisruptionFactory.convert(linxBreakends);
        LOGGER.info(" Loaded {} reportable disruptions from {}", reportableGeneDisruptions.size(), linxBreakendTsv);

                List<ReportableHomozygousDisruption> reportableHomozygousDisruptions =
                ReportableHomozygousDisruptionFactory.extractFromLinxDriverCatalogTsv(linxDriverCatalogTsv);
        LOGGER.info(" Loaded {} reportable homozygous disruptions from {}", reportableHomozygousDisruptions.size(), linxDriverCatalogTsv);

        return ImmutableLinxData.builder()
                .addAllFusions(linxFusions)
                .addAllGeneDisruptions(reportableGeneDisruptions)
                .addAllHomozygousDisruptions(reportableHomozygousDisruptions)
                .build();
    }
}
