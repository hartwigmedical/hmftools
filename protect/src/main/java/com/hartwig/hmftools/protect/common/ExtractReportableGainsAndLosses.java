package com.hartwig.hmftools.protect.common;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.drivercatalog.CNADrivers;
import com.hartwig.hmftools.common.drivercatalog.DriverCatalog;
import com.hartwig.hmftools.common.drivercatalog.panel.DriverGenePanel;
import com.hartwig.hmftools.common.purple.copynumber.CopyNumberInterpretation;
import com.hartwig.hmftools.common.purple.gene.GeneCopyNumber;

import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;

public class ExtractReportableGainsAndLosses {

    private ExtractReportableGainsAndLosses() {
    }

    @NotNull
    public static List<ReportableGainLoss> toReportableGainsAndLosses(@NotNull DriverGenePanel genePanel, @NotNull List<GeneCopyNumber> geneCopyNumbers, double ploidy) {
        CNADrivers copyNumberDriverModel = new CNADrivers(genePanel);
        List<DriverCatalog> drivers = Lists.newArrayList();
        drivers.addAll(copyNumberDriverModel.amplifications(ploidy, geneCopyNumbers));
        drivers.addAll(copyNumberDriverModel.deletions(geneCopyNumbers));

        List<ReportableGainLoss> reportableGainsAndLosses = Lists.newArrayList();
        for (DriverCatalog driver : drivers) {
            reportableGainsAndLosses.add(ImmutableReportableGainLoss.builder()
                    .chromosome(driver.chromosome())
                    .chromosomeBand(driver.chromosomeBand())
                    .gene(isCentromereOrTelomere(driver) ? Strings.EMPTY : driver.gene())
                    .interpretation(CopyNumberInterpretation.fromCNADriver(driver))
                    .copies(Math.round(Math.max(0, driver.minCopyNumber())))
                    .build());
        }

        return reportableGainsAndLosses;
    }

    private static boolean isCentromereOrTelomere(@NotNull DriverCatalog driver) {
        String band = driver.chromosomeBand().toLowerCase();
        return band.contains("centromere") || band.contains("telomere");
    }
}
