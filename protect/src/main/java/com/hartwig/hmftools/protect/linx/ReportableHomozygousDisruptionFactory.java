package com.hartwig.hmftools.protect.linx;

import java.io.IOException;
import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.drivercatalog.DriverCatalog;
import com.hartwig.hmftools.common.drivercatalog.DriverCatalogFile;
import com.hartwig.hmftools.common.drivercatalog.DriverType;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

public final class ReportableHomozygousDisruptionFactory {

    private static final Logger LOGGER = LogManager.getLogger(ReportableHomozygousDisruptionFactory.class);

    private ReportableHomozygousDisruptionFactory() {
    }

    @NotNull
    public static List<ReportableHomozygousDisruption> extractFromLinxDriversTsv(@NotNull String linxDriversTsv) throws IOException {
        List<DriverCatalog> linxDriversCatalog = DriverCatalogFile.read(linxDriversTsv);
        LOGGER.debug(" Loaded {} linx driver catalog records from {}", linxDriversCatalog.size(), linxDriversTsv);

        List<ReportableHomozygousDisruption> reportableHomozygousDisruptions = extractHomozygousDisruptions(linxDriversCatalog);
        LOGGER.debug("  Extracted {} homozygous disruptions from linx drivers", reportableHomozygousDisruptions.size());
        return reportableHomozygousDisruptions;
    }

    @NotNull
    private static List<ReportableHomozygousDisruption> extractHomozygousDisruptions(@NotNull List<DriverCatalog> linxDriversCatalog) {
        List<ReportableHomozygousDisruption> reportableHomozygousDisruptions = Lists.newArrayList();
        for (DriverCatalog driverCatalog : linxDriversCatalog) {
            if (driverCatalog.driver().equals(DriverType.HOM_DISRUPTION)) {
                reportableHomozygousDisruptions.add(ImmutableReportableHomozygousDisruption.builder()
                        .chromosome(driverCatalog.chromosome())
                        .chromosomeBand(driverCatalog.chromosomeBand())
                        .gene(driverCatalog.gene())
                        .build());
            }
        }

        return reportableHomozygousDisruptions;
    }
}
