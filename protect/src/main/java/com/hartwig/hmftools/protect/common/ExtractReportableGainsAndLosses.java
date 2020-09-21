package com.hartwig.hmftools.protect.common;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.drivercatalog.CNADrivers;
import com.hartwig.hmftools.common.drivercatalog.DriverCatalog;
import com.hartwig.hmftools.common.drivercatalog.DriverType;
import com.hartwig.hmftools.common.drivercatalog.panel.DriverGenePanel;
import com.hartwig.hmftools.common.purple.copynumber.CopyNumberInterpretation;
import com.hartwig.hmftools.common.purple.gene.GeneCopyNumber;

import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;

public class ExtractReportableGainsAndLosses {

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
                if (geneCopyNumber.germlineHet2HomRegions() > 0 || geneCopyNumber.germlineHomRegions() > 0) {
                    includeInReport = false;
                }
                includeInReport = true;
            }

            if (includeInReport) {
//                reportableGainsAndLosses.add(com.hartwig.hmftools.common.purple.copynumber.ImmutableReportableGainLoss.builder()
//                        .chromosome(driver.chromosome())
//                        .chromosomeBand(driver.chromosomeBand())
//                        .gene(isCentromereOrTelomere(driver) ? Strings.EMPTY : driver.gene())
//                        .interpretation(CopyNumberInterpretation.fromCNADriver(driver))
//                        .copies(Math.round(Math.max(0, driver.minCopyNumber())))
//                        .build());
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

    private static boolean isCentromereOrTelomere(@NotNull DriverCatalog driver) {
        String band = driver.chromosomeBand().toLowerCase();
        return band.contains("centromere") || band.contains("telomere");
    }
}
