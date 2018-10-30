package com.hartwig.hmftools.patientreporter.copynumber;

import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.actionability.ActionabilityAnalyzer;
import com.hartwig.hmftools.common.actionability.EvidenceItem;
import com.hartwig.hmftools.common.ecrf.projections.PatientTumorLocation;
import com.hartwig.hmftools.common.purple.copynumber.PurpleCopyNumber;
import com.hartwig.hmftools.common.purple.gene.GeneCopyNumber;
import com.hartwig.hmftools.common.purple.purity.FittedPurity;
import com.hartwig.hmftools.common.purple.purity.PurityContext;
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

        // KODU: Add gene copy numbers for which evidence has been found but which were not selected yet.
        for (GeneCopyNumber actionableGeneCopyNumber : evidencePerGeneCopyNumber.keySet()) {
            if (!evidencePerGeneCopyNumber.get(actionableGeneCopyNumber).isEmpty() && !reportableGeneCopyNumbers.contains(
                    actionableGeneCopyNumber)) {
                reportableGeneCopyNumbers.add(actionableGeneCopyNumber);
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
                .evidenceItems(toList(evidencePerGeneCopyNumber))
                .build();
    }

    @NotNull
    private static List<EvidenceItem> toList(@NotNull Map<GeneCopyNumber, List<EvidenceItem>> evidencePerCopyNumber) {
        List<EvidenceItem> evidenceItemList = Lists.newArrayList();
        for (List<EvidenceItem> items : evidencePerCopyNumber.values()) {
            evidenceItemList.addAll(items);
        }
        return evidenceItemList;
    }
}
