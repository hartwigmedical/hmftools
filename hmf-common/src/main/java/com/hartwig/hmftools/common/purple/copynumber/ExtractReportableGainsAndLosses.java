package com.hartwig.hmftools.common.purple.copynumber;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.drivercatalog.DriverCatalog;
import com.hartwig.hmftools.common.drivercatalog.DriverType;
import com.hartwig.hmftools.common.purple.gene.GeneCopyNumber;

import org.jetbrains.annotations.NotNull;

public final class ExtractReportableGainsAndLosses {

    private ExtractReportableGainsAndLosses() {
    }

    @NotNull
    public static List<ReportableGainLoss> toReportableGainsAndLosses(@NotNull List<DriverCatalog> driverCatalog,
            @NotNull List<GeneCopyNumber> geneCopyNumbers) {
        List<ReportableGainLoss> reportableGainsAndLosses = Lists.newArrayList();

        for (DriverCatalog driver : driverCatalog) {
            boolean includeInReport = false;
            if (driver.driver() == DriverType.AMP) {
                includeInReport = true;
            } else if (driver.driver() == DriverType.DEL) {
                GeneCopyNumber geneCopyNumber = findByGeneName(geneCopyNumbers, driver.gene());

                // Exclude DELs which are partially germline (see INC-34)
                includeInReport = geneCopyNumber.germlineHet2HomRegions() == 0 && geneCopyNumber.germlineHomRegions() == 0;
            }

            if (includeInReport) {
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

    @NotNull
    private static GeneCopyNumber findByGeneName(@NotNull List<GeneCopyNumber> geneCopyNumbers, @NotNull String gene) {
        for (GeneCopyNumber geneCopyNumber : geneCopyNumbers) {
            if (geneCopyNumber.gene().equals(gene)) {
                return geneCopyNumber;
            }
        }

        throw new IllegalStateException("Could not find gene across gene copy numbers: " + gene);
    }
}
