package com.hartwig.hmftools.common.linx;

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
    public static List<ReportableHomozygousDisruption> extractFromLinxDriverCatalogTsv(@NotNull String linxDriverCatalogTsv)
            throws IOException {
        List<DriverCatalog> linxDriversCatalog = DriverCatalogFile.read(linxDriverCatalogTsv);
        LOGGER.debug(" Loaded {} linx driver catalog records from {}", linxDriversCatalog.size(), linxDriverCatalogTsv);

        List<ReportableHomozygousDisruption> reportableHomozygousDisruptions = extractHomozygousDisruptions(linxDriversCatalog);
        LOGGER.debug("  Extracted {} homozygous disruptions from linx drivers", reportableHomozygousDisruptions.size());
        return reportableHomozygousDisruptions;
    }

    @NotNull
    private static List<ReportableHomozygousDisruption> extractHomozygousDisruptions(@NotNull List<DriverCatalog> driverCatalog) {
        List<ReportableHomozygousDisruption> homozygousDisruptions = Lists.newArrayList();

        for (DriverCatalog driver : driverCatalog) {
            if (driver.driver() == DriverType.HOM_DUP_DISRUPTION || driver.driver() == DriverType.HOM_DEL_DISRUPTION) {
                homozygousDisruptions.add(create(driver));
            }
        }

        return homozygousDisruptions;
    }

    @NotNull
    private static ReportableHomozygousDisruption create(@NotNull DriverCatalog driverCatalog) {
        return ImmutableReportableHomozygousDisruption.builder()
                .chromosome(driverCatalog.chromosome())
                .chromosomeBand(driverCatalog.chromosomeBand())
                .gene(driverCatalog.gene())
                .transcript(driverCatalog.transcript())
                .isCanonical(driverCatalog.isCanonical())
                .build();
    }
}
