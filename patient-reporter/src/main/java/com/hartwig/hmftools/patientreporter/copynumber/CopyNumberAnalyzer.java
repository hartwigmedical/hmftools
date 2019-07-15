package com.hartwig.hmftools.patientreporter.copynumber;

import java.util.Collections;
import java.util.List;
import java.util.Map;

import com.hartwig.hmftools.common.actionability.ActionabilityAnalyzer;
import com.hartwig.hmftools.common.actionability.EvidenceItem;
import com.hartwig.hmftools.common.ecrf.projections.PatientTumorLocation;
import com.hartwig.hmftools.common.purple.gene.GeneCopyNumber;
import com.hartwig.hmftools.common.purple.purity.FittedPurity;
import com.hartwig.hmftools.common.purple.purity.FittedPurityStatus;
import com.hartwig.hmftools.common.purple.purity.PurityContext;
import com.hartwig.hmftools.patientreporter.actionability.ReportableEvidenceItemFactory;
import com.hartwig.hmftools.patientreporter.driver.DriverGeneView;

import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public final class CopyNumberAnalyzer {

    private CopyNumberAnalyzer() {
    }

    @NotNull
    public static CopyNumberAnalysis run(@NotNull PurityContext purityContext, @NotNull List<GeneCopyNumber> exomeGeneCopyNumbers,
            @NotNull DriverGeneView driverGeneView, @NotNull ActionabilityAnalyzer actionabilityAnalyzer,
            @Nullable PatientTumorLocation patientTumorLocation) {
        FittedPurity bestFit = purityContext.bestFit();

        List<GeneCopyNumber> reportableGeneCopyNumbers =
                ReportingCopyNumberFilters.filterForReporting(exomeGeneCopyNumbers, driverGeneView, purityContext.gender(), bestFit.ploidy());

        String primaryTumorLocation = patientTumorLocation != null ? patientTumorLocation.primaryTumorLocation() : null;
        Map<GeneCopyNumber, List<EvidenceItem>> evidencePerGeneCopyNumber =
                actionabilityAnalyzer.evidenceForCopyNumbers(exomeGeneCopyNumbers, primaryTumorLocation, bestFit.ploidy());

        List<EvidenceItem> filteredEvidence = ReportableEvidenceItemFactory.reportableFlatList(evidencePerGeneCopyNumber);

        // Add gene copy numbers for which filtered evidence has been found but which were not selected yet.
        for (Map.Entry<GeneCopyNumber, List<EvidenceItem>> entry : evidencePerGeneCopyNumber.entrySet()) {
            GeneCopyNumber geneCopyNumber = entry.getKey();
            if (!Collections.disjoint(entry.getValue(), filteredEvidence) && !reportableGeneCopyNumbers.contains(geneCopyNumber)) {
                reportableGeneCopyNumbers.add(geneCopyNumber);
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
}
