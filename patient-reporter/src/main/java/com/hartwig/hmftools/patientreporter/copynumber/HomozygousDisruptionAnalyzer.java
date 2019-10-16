package com.hartwig.hmftools.patientreporter.copynumber;

import java.io.IOException;
import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.drivercatalog.DriverCatalog;
import com.hartwig.hmftools.common.drivercatalog.DriverCatalogFile;
import com.hartwig.hmftools.common.drivercatalog.DriverType;
import com.hartwig.hmftools.common.drivercatalog.LikelihoodMethod;
import com.hartwig.hmftools.common.variant.structural.linx.ImmutableLinxDriver;
import com.hartwig.hmftools.common.variant.structural.linx.LinxDriver;
import com.hartwig.hmftools.common.variant.structural.linx.LinxDriverFile;
import com.hartwig.hmftools.patientreporter.structural.ImmutableReportableDriverCatalog;
import com.hartwig.hmftools.patientreporter.structural.ReportableDriverCatalog;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

public final class HomozygousDisruptionAnalyzer {

    private static final Logger LOGGER = LogManager.getLogger(HomozygousDisruptionAnalyzer.class);

    private HomozygousDisruptionAnalyzer() {

    }

    @NotNull
    public static List<ReportableDriverCatalog> interpetDriverCatalog(@NotNull String linxDriversCatalogTsv) throws IOException {
        List<DriverCatalog> allDriversCatalog = DriverCatalogFile.read(linxDriversCatalogTsv);
        LOGGER.info("Loaded {} driver catalog records", allDriversCatalog.size());

        List<ReportableDriverCatalog> reportableDriverCatalogs = extractHomozygousDisruptions(allDriversCatalog);
        LOGGER.info("Loaded {} homozyguous deletion disruptions records", reportableDriverCatalogs.size());
        return reportableDriverCatalogs;
    }

    @NotNull
    public static List<ReportableDriverCatalog> extractHomozygousDisruptions(List<DriverCatalog> allDriversCatalog) {
        List<ReportableDriverCatalog> reportableDriverCatalogs = Lists.newArrayList();
        for (DriverCatalog driverCatalog : allDriversCatalog) {
            if (driverCatalog.driver().equals(DriverType.HOM_DISRUPTION) && driverCatalog.likelihoodMethod().equals(LikelihoodMethod.DEL)) {
                reportableDriverCatalogs.add(ImmutableReportableDriverCatalog.builder()
                        .chromosome(driverCatalog.chromosome())
                        .chromosomeBand(driverCatalog.chromosomeBand())
                        .gene(driverCatalog.gene())
                        .driver(driverCatalog.driver().toString().toLowerCase())
                        .build());
            }
        }
        return reportableDriverCatalogs;
    }
}
