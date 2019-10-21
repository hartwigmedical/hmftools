package com.hartwig.hmftools.patientreporter.copynumber;

import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.Set;

import com.google.common.collect.Lists;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.actionability.ActionabilityAnalyzer;
import com.hartwig.hmftools.common.actionability.EvidenceItem;
import com.hartwig.hmftools.common.drivercatalog.CNADrivers;
import com.hartwig.hmftools.common.drivercatalog.DriverCatalog;
import com.hartwig.hmftools.common.ecrf.projections.PatientTumorLocation;
import com.hartwig.hmftools.common.purple.gene.GeneCopyNumber;
import com.hartwig.hmftools.common.purple.purity.FittedPurity;
import com.hartwig.hmftools.common.purple.purity.FittedPurityStatus;
import com.hartwig.hmftools.common.purple.purity.PurityContext;
import com.hartwig.hmftools.patientreporter.actionability.ReportableEvidenceItemFactory;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public final class CopyNumberAnalyzer {

    private static final Logger LOGGER = LogManager.getLogger(CopyNumberAnalyzer.class);

    private CopyNumberAnalyzer() {
    }

    @NotNull
    public static CopyNumberAnalysis run(@NotNull PurityContext purityContext, @NotNull List<GeneCopyNumber> exomeGeneCopyNumbers,
            @NotNull ActionabilityAnalyzer actionabilityAnalyzer, @Nullable PatientTumorLocation patientTumorLocation) {
        FittedPurity bestFit = purityContext.bestFit();

        CNADrivers copyNumberDriverModel = new CNADrivers();
        List<DriverCatalog> drivers = Lists.newArrayList();
        drivers.addAll(copyNumberDriverModel.amplifications(bestFit.ploidy(), exomeGeneCopyNumbers));
        drivers.addAll(copyNumberDriverModel.deletions(exomeGeneCopyNumbers));

        String primaryTumorLocation = patientTumorLocation != null ? patientTumorLocation.primaryTumorLocation() : null;
        Map<GeneCopyNumber, List<EvidenceItem>> evidencePerGeneCopyNumber =
                actionabilityAnalyzer.evidenceForCopyNumbers(exomeGeneCopyNumbers, primaryTumorLocation, bestFit.ploidy());

        List<EvidenceItem> filteredEvidence = ReportableEvidenceItemFactory.reportableFlatList(evidencePerGeneCopyNumber);

        List<ReportableGainLoss> reportableGainsAndLosses = toReportableGainsAndLosses(drivers);

        // Check that all copy numbers with evidence are reported (since they are in the driver catalog).
        Set<String> reportableGenes = Sets.newHashSet();
        for (ReportableGainLoss gainLoss : reportableGainsAndLosses) {
            reportableGenes.add(gainLoss.gene());
        }

        for (Map.Entry<GeneCopyNumber, List<EvidenceItem>> entry : evidencePerGeneCopyNumber.entrySet()) {
            GeneCopyNumber geneCopyNumber = entry.getKey();
            if (!Collections.disjoint(entry.getValue(), filteredEvidence) && !reportableGenes.contains(geneCopyNumber.gene())) {
                LOGGER.warn("Copy number with evidence not reported: {}!", geneCopyNumber.gene());
            }
        }

        return ImmutableCopyNumberAnalysis.builder()
                .hasReliablePurityFit(purityContext.status() != FittedPurityStatus.NO_TUMOR)
                .purity(bestFit.purity())
                .ploidy(bestFit.ploidy())
                .exomeGeneCopyNumbers(exomeGeneCopyNumbers)
                .reportableGainsAndLosses(reportableGainsAndLosses)
                .evidenceItems(filteredEvidence)
                .build();
    }

    @NotNull
    private static List<ReportableGainLoss> toReportableGainsAndLosses(@NotNull List<DriverCatalog> drivers) {
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
