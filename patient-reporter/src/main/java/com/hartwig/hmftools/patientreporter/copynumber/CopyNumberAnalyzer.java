package com.hartwig.hmftools.patientreporter.copynumber;

import java.util.Collections;
import java.util.List;
import java.util.Map;

import com.hartwig.hmftools.common.actionability.ActionabilityAnalyzer;
import com.hartwig.hmftools.common.actionability.EvidenceItem;
import com.hartwig.hmftools.common.ecrf.projections.PatientTumorLocation;
import com.hartwig.hmftools.common.purple.copynumber.PurpleCopyNumber;
import com.hartwig.hmftools.common.purple.gene.GeneCopyNumber;
import com.hartwig.hmftools.common.purple.purity.FittedPurity;
import com.hartwig.hmftools.common.purple.purity.PurityContext;
import com.hartwig.hmftools.patientreporter.actionability.ReportableEvidenceItemFactory;
import com.hartwig.hmftools.patientreporter.genepanel.GeneModel;

import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public final class CopyNumberAnalyzer {

    private CopyNumberAnalyzer() {
    }

    @NotNull
    public static CopyNumberAnalysis run(@NotNull PurityContext purityContext, @NotNull List<PurpleCopyNumber> purpleCopyNumbers,
            @NotNull List<GeneCopyNumber> exomeGeneCopyNumbers, @NotNull GeneModel geneModel,
            @NotNull ActionabilityAnalyzer actionabilityAnalyzer, @Nullable PatientTumorLocation patientTumorLocation) {
        FittedPurity bestFit = purityContext.bestFit();

        List<GeneCopyNumber> significantCopyNumbers =
                ReportingCopyNumberFilters.filterForSignificance(exomeGeneCopyNumbers, bestFit.ploidy());

        List<GeneCopyNumber> reportableGeneCopyNumbers = ReportingCopyNumberFilters.filterForReporting(significantCopyNumbers, geneModel);

        String primaryTumorLocation = patientTumorLocation != null ? patientTumorLocation.primaryTumorLocation() : null;
        Map<GeneCopyNumber, List<EvidenceItem>> evidencePerGeneCopyNumber =
                actionabilityAnalyzer.evidenceForCopyNumbers(significantCopyNumbers, primaryTumorLocation);

        List<EvidenceItem> filteredEvidence = ReportableEvidenceItemFactory.reportableFlatList(evidencePerGeneCopyNumber);

        // KODU: Add gene copy numbers for which filtered evidence has been found but which were not selected yet.
        for (Map.Entry<GeneCopyNumber, List<EvidenceItem>> entry : evidencePerGeneCopyNumber.entrySet()) {
            GeneCopyNumber geneCopyNumber = entry.getKey();
            if (!Collections.disjoint(entry.getValue(), filteredEvidence) && !reportableGeneCopyNumbers.contains(geneCopyNumber)) {
                reportableGeneCopyNumbers.add(geneCopyNumber);
            }
        }

        return ImmutableCopyNumberAnalysis.builder()
                .gender(purityContext.gender())
                .status(purityContext.status())
                .fittedPurity(bestFit)
                .fittedScorePurity(purityContext.score())
                .copyNumbers(purpleCopyNumbers)
                .exomeGeneCopyNumbers(exomeGeneCopyNumbers)
                .reportableGeneCopyNumbers(reportableGeneCopyNumbers)
                .evidenceItems(filteredEvidence)
                .build();
    }
}
