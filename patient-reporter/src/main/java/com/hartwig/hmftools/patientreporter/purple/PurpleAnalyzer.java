package com.hartwig.hmftools.patientreporter.purple;

import java.util.List;
import java.util.Map;

import com.hartwig.hmftools.common.actionability.ActionabilityAnalyzer;
import com.hartwig.hmftools.common.actionability.EvidenceItem;
import com.hartwig.hmftools.common.drivercatalog.DriverCatalog;
import com.hartwig.hmftools.common.ecrf.projections.PatientTumorLocation;
import com.hartwig.hmftools.common.purple.CheckPurpleQuality;
import com.hartwig.hmftools.common.purple.copynumber.ExtractReportableGainsAndLosses;
import com.hartwig.hmftools.common.purple.copynumber.ReportableGainLoss;
import com.hartwig.hmftools.common.purple.gene.GeneCopyNumber;
import com.hartwig.hmftools.common.purple.purity.FittedPurity;
import com.hartwig.hmftools.common.purple.purity.PurityContext;
import com.hartwig.hmftools.common.purple.qc.PurpleQC;
import com.hartwig.hmftools.patientreporter.actionability.ReportableEvidenceItemFactory;

import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public final class PurpleAnalyzer {

    private PurpleAnalyzer() {
    }

    @NotNull
    public static PurpleAnalysis run(@NotNull PurityContext purityContext, @NotNull PurpleQC purpleQC,
            @NotNull List<GeneCopyNumber> exomeGeneCopyNumbers, @NotNull ActionabilityAnalyzer actionabilityAnalyzer,
            @Nullable PatientTumorLocation patientTumorLocation, @NotNull List<DriverCatalog> driverCatalog) {
        FittedPurity bestFit = purityContext.bestFit();
        List<ReportableGainLoss> reportableGainsAndLosses =
                ExtractReportableGainsAndLosses.toReportableGainsAndLosses(driverCatalog, exomeGeneCopyNumbers);
        String primaryTumorLocation = patientTumorLocation != null ? patientTumorLocation.primaryTumorLocation() : null;
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
                .hasReliableQuality(CheckPurpleQuality.checkHasReliableQuality(purpleQC))
                .ploidy(bestFit.ploidy())
                .exomeGeneCopyNumbers(exomeGeneCopyNumbers)
                .reportableGainsAndLosses(reportableGainsAndLosses)
                .purpleSignatures(purpleSignatures)
                .evidenceItems(filteredEvidenceItems)
                .build();
    }
}
