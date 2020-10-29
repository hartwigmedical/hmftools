package com.hartwig.hmftools.common.purple.copynumber;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.drivercatalog.DriverCatalog;
import com.hartwig.hmftools.common.drivercatalog.DriverType;

import org.jetbrains.annotations.NotNull;

public final class ExtractReportableGainsAndLosses {

    private ExtractReportableGainsAndLosses() {
    }

    @NotNull
    public static List<ReportableGainLoss> toReportableGainsAndLosses(@NotNull List<DriverCatalog> driverCatalog) {
        List<ReportableGainLoss> reportableGainsAndLosses = Lists.newArrayList();

        for (DriverCatalog driver : driverCatalog) {
            if (driver.driver() == DriverType.AMP || driver.driver() == DriverType.DEL) {
                reportableGainsAndLosses.add(ImmutableReportableGainLoss.builder()
                        .chromosome(driver.chromosome())
                        .chromosomeBand(driver.chromosomeBand())
                        .gene(driver.gene())
                        .interpretation(CopyNumberInterpretation.fromCNADriver(driver))
                        .copies(Math.round(Math.max(0, driver.minCopyNumber())))
                        .build());
            }
        }

        return reportableGainsAndLosses;
    }
}
