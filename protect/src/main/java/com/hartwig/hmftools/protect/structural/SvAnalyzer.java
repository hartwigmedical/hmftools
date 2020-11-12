package com.hartwig.hmftools.protect.structural;

import java.util.List;
import java.util.Map;

import com.hartwig.hmftools.common.actionability.ActionabilityAnalyzer;
import com.hartwig.hmftools.common.actionability.EvidenceItem;
import com.hartwig.hmftools.common.clinical.PatientPrimaryTumor;
import com.hartwig.hmftools.common.variant.structural.linx.LinxBreakend;
import com.hartwig.hmftools.common.variant.structural.linx.LinxFusion;
import com.hartwig.hmftools.protect.actionability.ReportableEvidenceItemFactory;
import com.hartwig.hmftools.protect.linx.LinxData;

import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public final class SvAnalyzer {

    private SvAnalyzer() {
    }

    @NotNull
    public static SvAnalysis run(@NotNull LinxData linxData, @NotNull ActionabilityAnalyzer actionabilityAnalyzer,
            @Nullable PatientPrimaryTumor patientPrimaryTumor) {
        String primaryTumorLocation = patientPrimaryTumor != null ? patientPrimaryTumor.primaryTumorLocation() : null;
        Map<LinxFusion, List<EvidenceItem>> evidencePerFusion =
                actionabilityAnalyzer.evidenceForFusions(linxData.fusions(), primaryTumorLocation);

        List<EvidenceItem> filteredEvidence = ReportableEvidenceItemFactory.toReportableFlatList(evidencePerFusion);

        return ImmutableSvAnalysis.builder()
                .reportableFusions(linxData.fusions())
                .reportableDisruptions(linxData.geneDisruptions())
                .evidenceItems(filteredEvidence)
                .build();
    }

    @NotNull
    public static SvAnalysis run(@NotNull List<LinxFusion> fusions, @NotNull List<LinxBreakend> disruptions,
            @NotNull ActionabilityAnalyzer actionabilityAnalyzer, @Nullable PatientPrimaryTumor patientPrimaryTumor) {
        List<ReportableGeneDisruption> reportableGeneDisruptions = ReportableGeneDisruptionFactory.convert(disruptions);

        String primaryTumorLocation = patientPrimaryTumor != null ? patientPrimaryTumor.primaryTumorLocation() : null;

        Map<LinxFusion, List<EvidenceItem>> evidencePerFusion = actionabilityAnalyzer.evidenceForFusions(fusions, primaryTumorLocation);

        List<EvidenceItem> filteredEvidence = ReportableEvidenceItemFactory.toReportableFlatList(evidencePerFusion);

        return ImmutableSvAnalysis.builder()
                .reportableFusions(fusions)
                .reportableDisruptions(reportableGeneDisruptions)
                .evidenceItems(filteredEvidence)
                .build();
    }
}
