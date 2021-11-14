package com.hartwig.hmftools.common.linx;

import java.io.File;
import java.io.IOException;
import java.util.List;
import java.util.stream.Collectors;

import com.hartwig.hmftools.common.sv.linx.LinxBreakend;
import com.hartwig.hmftools.common.sv.linx.LinxDriver;
import com.hartwig.hmftools.common.sv.linx.LinxFusion;

import org.apache.commons.compress.utils.Lists;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public final class LinxDataLoader {

    private static final Logger LOGGER = LogManager.getLogger(LinxDataLoader.class);

    private LinxDataLoader() {
    }

    @NotNull
    public static LinxData load(@NotNull String linxFusionTsv, @NotNull String linxBreakendTsv, @NotNull String linxDriverCatalogTsv)
            throws IOException {
       return load(linxFusionTsv, linxBreakendTsv, linxDriverCatalogTsv, null);
    }

    @NotNull
    public static LinxData load(@NotNull String linxFusionTsv, @NotNull String linxBreakendTsv, @NotNull String linxDriverCatalogTsv,
            @Nullable String linxDriverTsv) throws IOException {
        LOGGER.info("Loading LINX data from {}", new File(linxFusionTsv).getParent());
        List<LinxFusion> fusions = LinxFusion.read(linxFusionTsv);

        List<LinxFusion> reportableFusions = Lists.newArrayList();
        List<LinxFusion> unreportedFusions = Lists.newArrayList();
        for (LinxFusion fusion : fusions) {
            if (fusion.reported()) {
                reportableFusions.add(fusion);
            } else {
                unreportedFusions.add(fusion);
            }
        }
        LOGGER.info(" Loaded {} fusions (of which {} are reportable) from {}", fusions.size(), reportableFusions.size(), linxFusionTsv);

        List<LinxBreakend> linxBreakends =
                LinxBreakend.read(linxBreakendTsv).stream().filter(LinxBreakend::reportedDisruption).collect(Collectors.toList());
        List<ReportableGeneDisruption> geneDisruptions = ReportableGeneDisruptionFactory.convert(linxBreakends);
        LOGGER.debug(" Generated {} reportable disruptions based on {} breakends", geneDisruptions.size(), linxBreakends.size());
        LOGGER.info(" Loaded {} reportable disruptions from {}", geneDisruptions.size(), linxBreakendTsv);

        List<ReportableHomozygousDisruption> homozygousDisruptions =
                ReportableHomozygousDisruptionFactory.extractFromLinxDriverCatalogTsv(linxDriverCatalogTsv);
        LOGGER.info(" Loaded {} reportable homozygous disruptions from {}", homozygousDisruptions.size(), linxDriverCatalogTsv);

        List<LinxDriver> drivers = Lists.newArrayList();
        if (linxDriverTsv != null) {
            drivers = LinxDriver.read(linxDriverTsv);
            LOGGER.info(" Loaded {} drivers from {}", drivers.size(), linxDriverTsv);
        }

        return ImmutableLinxData.builder()
                .reportableFusions(reportableFusions)
                .unreportedFusions(unreportedFusions)
                .geneDisruptions(geneDisruptions)
                .homozygousDisruptions(homozygousDisruptions)
                .drivers(drivers)
                .build();
    }
}
