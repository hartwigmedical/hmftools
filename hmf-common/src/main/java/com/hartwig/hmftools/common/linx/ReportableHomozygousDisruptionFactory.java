package com.hartwig.hmftools.common.linx;

import java.io.IOException;
import java.util.List;
import java.util.Map;
import java.util.Set;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.drivercatalog.DriverCatalog;
import com.hartwig.hmftools.common.drivercatalog.DriverCatalogFile;
import com.hartwig.hmftools.common.drivercatalog.DriverCatalogKey;
import com.hartwig.hmftools.common.drivercatalog.DriverCatalogMap;
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
    private static List<ReportableHomozygousDisruption> extractHomozygousDisruptions(@NotNull List<DriverCatalog> linxDriversCatalog) {
        List<ReportableHomozygousDisruption> homozygousDisruptions = Lists.newArrayList();
        Set<DriverCatalogKey> keys = DriverCatalogKey.buildUniqueKeysSet(linxDriversCatalog);
        Map<DriverCatalogKey, DriverCatalog> geneDriverMap = DriverCatalogMap.toDriverMap(linxDriversCatalog);

        for (DriverCatalogKey key : keys) {
            DriverCatalog geneDriver = geneDriverMap.get(key);
            if (geneDriver == null) {
                throw new IllegalStateException(
                        "Could not find driver entry for homozygous disruption on gene '" + geneDriver.gene() + "'");
            }
            if (geneDriver.driver() == DriverType.HOM_DUP_DISRUPTION || geneDriver.driver() == DriverType.HOM_DEL_DISRUPTION) {
                homozygousDisruptions.add(create(geneDriver));
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
