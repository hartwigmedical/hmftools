package com.hartwig.hmftools.patientreporter.copynumber;

import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.stream.Collectors;

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

        CNADrivers copyNumberDrivers = new CNADrivers();
        List<DriverCatalog> drivers = Lists.newArrayList();
        drivers.addAll(copyNumberDrivers.amplifications(bestFit.ploidy(), exomeGeneCopyNumbers));
        drivers.addAll(copyNumberDrivers.deletions(exomeGeneCopyNumbers));

        String primaryTumorLocation = patientTumorLocation != null ? patientTumorLocation.primaryTumorLocation() : null;
        Map<GeneCopyNumber, List<EvidenceItem>> evidencePerGeneCopyNumber =
                actionabilityAnalyzer.evidenceForCopyNumbers(exomeGeneCopyNumbers, primaryTumorLocation, bestFit.ploidy());

        List<EvidenceItem> filteredEvidence = ReportableEvidenceItemFactory.reportableFlatList(evidencePerGeneCopyNumber);

        List<GeneCopyNumber> reportableGeneCopyNumbers = extractReportableCopyNumbers(exomeGeneCopyNumbers, drivers);

        // Check that all copy numbers with evidence are reported (since they are in the driver catalog).
        for (Map.Entry<GeneCopyNumber, List<EvidenceItem>> entry : evidencePerGeneCopyNumber.entrySet()) {
            GeneCopyNumber geneCopyNumber = entry.getKey();
            if (!Collections.disjoint(entry.getValue(), filteredEvidence) && !reportableGeneCopyNumbers.contains(geneCopyNumber)) {
                LOGGER.warn("Copy number with evidence not reported: {}!", geneCopyNumber.gene());
            }
        }

        return ImmutableCopyNumberAnalysis.builder()
                .hasReliablePurityFit(purityContext.status() != FittedPurityStatus.NO_TUMOR)
                .purity(bestFit.purity())
                .ploidy(bestFit.ploidy())
                .gender(purityContext.gender())
                .exomeGeneCopyNumbers(exomeGeneCopyNumbers)
                .reportableGeneCopyNumbers(reportableGeneCopyNumbers)
                .evidenceItems(filteredEvidence)
                .build();
    }

    @NotNull
    private static List<GeneCopyNumber> extractReportableCopyNumbers(@NotNull List<GeneCopyNumber> exomeGeneCopyNumbers,
            @NotNull List<DriverCatalog> drivers) {
        Set<String> driverGenes = Sets.newHashSet();
        for (DriverCatalog driver : drivers) {
            driverGenes.add(driver.gene());
        }

        return exomeGeneCopyNumbers.stream().filter(copyNumber -> driverGenes.contains(copyNumber.gene())).collect(Collectors.toList());
    }
}
