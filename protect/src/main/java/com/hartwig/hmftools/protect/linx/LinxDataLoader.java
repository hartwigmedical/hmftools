package com.hartwig.hmftools.protect.linx;

import java.io.File;
import java.io.IOException;
import java.util.List;
import java.util.stream.Collectors;

import com.hartwig.hmftools.common.variant.structural.linx.LinxBreakend;
import com.hartwig.hmftools.common.variant.structural.linx.LinxFusion;
import com.hartwig.hmftools.protect.ProtectConfig;

import org.apache.commons.compress.utils.Lists;
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
    public static LinxData load(@NotNull String linxFusionTsv, @NotNull String linxBreakendTsv, @NotNull String linxDriverCatalogTsv)
            throws IOException {
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

        return ImmutableLinxData.builder()
                .reportableFusions(reportableFusions)
                .unreportedFusions(unreportedFusions)
                .geneDisruptions(geneDisruptions)
                .homozygousDisruptions(homozygousDisruptions)
                .build();
    }
}
