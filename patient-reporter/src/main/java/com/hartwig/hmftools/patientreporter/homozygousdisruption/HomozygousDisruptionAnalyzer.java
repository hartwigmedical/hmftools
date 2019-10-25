package com.hartwig.hmftools.patientreporter.homozygousdisruption;

import java.io.IOException;
import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.drivercatalog.DriverCatalog;
import com.hartwig.hmftools.common.drivercatalog.DriverCatalogFile;
import com.hartwig.hmftools.common.drivercatalog.DriverType;
import com.hartwig.hmftools.common.drivercatalog.LikelihoodMethod;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

public final class HomozygousDisruptionAnalyzer {

    private static final Logger LOGGER = LogManager.getLogger(HomozygousDisruptionAnalyzer.class);

    private HomozygousDisruptionAnalyzer() {

    }

    @NotNull
    public static List<ReportableHomozygousDisruption> interpetDriverCatalog(@NotNull String linxDriversCatalogTsv) throws IOException {
        List<DriverCatalog> allDriversCatalog = DriverCatalogFile.read(linxDriversCatalogTsv);
        LOGGER.info("Loaded {} driver catalog records", allDriversCatalog.size());

        List<ReportableHomozygousDisruption> reportableHomozygousDisruptions = extractHomozygousDisruptions(allDriversCatalog);
        LOGGER.info("Loaded {} homozyguous deletion disruptions records", reportableHomozygousDisruptions.size());
        return reportableHomozygousDisruptions;
    }

    @NotNull
    public static List<ReportableHomozygousDisruption> extractHomozygousDisruptions(List<DriverCatalog> allDriversCatalog) {
        List<ReportableHomozygousDisruption> reportableHomozygousDisruptions = Lists.newArrayList();
        for (DriverCatalog driverCatalog : allDriversCatalog) {
            if (driverCatalog.driver().equals(DriverType.HOM_DISRUPTION) && driverCatalog.likelihoodMethod().equals(LikelihoodMethod.DEL)) {
                reportableHomozygousDisruptions.add(ImmutableReportableHomozygousDisruption.builder()
                        .chromosome(driverCatalog.chromosome())
                        .chromosomeBand(driverCatalog.chromosomeBand())
                        .gene(driverCatalog.gene())
                        .driver(driverCatalog.driver().toString().toLowerCase())
                        .build());
            }
        }
        return reportableHomozygousDisruptions;
    }
}
