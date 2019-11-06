package com.hartwig.hmftools.patientreporter.copynumber;

import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.Set;

import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.actionability.ActionabilityAnalyzer;
import com.hartwig.hmftools.common.actionability.EvidenceItem;
import com.hartwig.hmftools.common.reportablegenomicalterations.ReportableGainLoss;
import com.hartwig.hmftools.common.reportablegenomicalterations.extractReportableGainsAndLosses;
import com.hartwig.hmftools.common.ecrf.projections.PatientTumorLocation;
import com.hartwig.hmftools.common.purple.gene.GeneCopyNumber;
import com.hartwig.hmftools.common.purple.purity.FittedPurity;
import com.hartwig.hmftools.common.purple.purity.FittedPurityStatus;
import com.hartwig.hmftools.common.purple.purity.PurityContext;
import com.hartwig.hmftools.common.purple.qc.PurpleQC;
import com.hartwig.hmftools.common.purple.qc.PurpleQCStatus;
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
    public static CopyNumberAnalysis run(@NotNull PurityContext purityContext, @NotNull PurpleQC purpleQC,
            @NotNull List<GeneCopyNumber> exomeGeneCopyNumbers, @NotNull ActionabilityAnalyzer actionabilityAnalyzer,
            @Nullable PatientTumorLocation patientTumorLocation) {
        FittedPurity bestFit = purityContext.bestFit();
        List<ReportableGainLoss> reportableGainsAndLosses =
                extractReportableGainsAndLosses.toReportableGainsAndLosses(exomeGeneCopyNumbers, bestFit.ploidy());
        String primaryTumorLocation = patientTumorLocation != null ? patientTumorLocation.primaryTumorLocation() : null;
        Map<GeneCopyNumber, List<EvidenceItem>> evidencePerGeneCopyNumber =
                actionabilityAnalyzer.evidenceForCopyNumbers(exomeGeneCopyNumbers, primaryTumorLocation, bestFit.ploidy());

        List<EvidenceItem> filteredEvidence = ReportableEvidenceItemFactory.toReportableFlatList(evidencePerGeneCopyNumber);

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
                .purity(bestFit.purity())
                .hasReliablePurity(purityContext.status() != FittedPurityStatus.NO_TUMOR)
                .hasReliableQuality(purpleQC.status() == PurpleQCStatus.PASS)
                .ploidy(bestFit.ploidy())
                .exomeGeneCopyNumbers(exomeGeneCopyNumbers)
                .reportableGainsAndLosses(reportableGainsAndLosses)
                .evidenceItems(filteredEvidence)
                .build();
    }
}
