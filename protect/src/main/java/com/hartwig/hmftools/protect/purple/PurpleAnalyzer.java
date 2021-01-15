package com.hartwig.hmftools.protect.purple;

import java.util.List;
import java.util.Map;

import com.hartwig.hmftools.common.actionability.ActionabilityAnalyzer;
import com.hartwig.hmftools.common.actionability.EvidenceItem;
import com.hartwig.hmftools.common.clinical.PatientPrimaryTumor;
import com.hartwig.hmftools.common.drivercatalog.DriverCatalog;
import com.hartwig.hmftools.common.purple.CheckPurpleQuality;
import com.hartwig.hmftools.common.purple.copynumber.ExtractReportableGainsAndLosses;
import com.hartwig.hmftools.common.purple.copynumber.ReportableGainLoss;
import com.hartwig.hmftools.common.purple.purity.FittedPurity;
import com.hartwig.hmftools.common.purple.purity.PurityContext;
import com.hartwig.hmftools.common.purple.qc.PurpleQC;
import com.hartwig.hmftools.protect.actionability.ReportableEvidenceItemFactory;

import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public final class PurpleAnalyzer {

    private PurpleAnalyzer() {
    }

    @NotNull
    public static List<EvidenceItem> run(@NotNull PurpleData purpleData,
            @NotNull ActionabilityAnalyzer actionabilityAnalyzer, @Nullable PatientPrimaryTumor patientPrimaryTumor) {
        String primaryTumorLocation = patientPrimaryTumor != null ? patientPrimaryTumor.location() : null;
        Map<ReportableGainLoss, List<EvidenceItem>> evidencePerGeneCopyNumber =
                actionabilityAnalyzer.evidenceForCopyNumbers(purpleData.copyNumberAlterations(), primaryTumorLocation);

        return ReportableEvidenceItemFactory.toReportableFlatList(evidencePerGeneCopyNumber);
    }

    @NotNull
    public static PurpleAnalysis run(@NotNull PurityContext purityContext, @NotNull PurpleQC purpleQC,
            @NotNull ActionabilityAnalyzer actionabilityAnalyzer, @Nullable PatientPrimaryTumor patientPrimaryTumor,
            @NotNull List<DriverCatalog> purpleDriverCatalog) {
        FittedPurity bestFit = purityContext.bestFit();
        List<ReportableGainLoss> reportableGainsAndLosses = ExtractReportableGainsAndLosses.toReportableGainsAndLosses(purpleDriverCatalog);
        String primaryTumorLocation = patientPrimaryTumor != null ? patientPrimaryTumor.location() : null;
        Map<ReportableGainLoss, List<EvidenceItem>> evidencePerGeneCopyNumber =
                actionabilityAnalyzer.evidenceForCopyNumbers(reportableGainsAndLosses, primaryTumorLocation);

        List<EvidenceItem> filteredEvidenceItems = ReportableEvidenceItemFactory.toReportableFlatList(evidencePerGeneCopyNumber);

        PurpleSignatures purpleSignatures = ImmutablePurpleSignatures.builder()
                .microsatelliteIndelsPerMb(purityContext.microsatelliteIndelsPerMb())
                .microsatelliteStatus(purityContext.microsatelliteStatus())
                .tumorMutationalBurdenPerMb(purityContext.tumorMutationalBurdenPerMb())
                .tumorMutationalLoad(purityContext.tumorMutationalLoad())
                .tumorMutationalLoadStatus(purityContext.tumorMutationalLoadStatus())
                .build();

        return ImmutablePurpleAnalysis.builder()
                .purity(bestFit.purity())
                .hasReliablePurity(CheckPurpleQuality.checkHasReliablePurity(purityContext))
                .hasReliableQuality(purpleQC.pass())
                .ploidy(bestFit.ploidy())
                .reportableGainsAndLosses(reportableGainsAndLosses)
                .purpleSignatures(purpleSignatures)
                .evidenceItems(filteredEvidenceItems)
                .build();
    }
}
