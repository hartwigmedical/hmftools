package com.hartwig.hmftools.patientreporter.structural;

import java.util.List;
import java.util.Map;

import com.hartwig.hmftools.common.actionability.ActionabilityAnalyzer;
import com.hartwig.hmftools.common.actionability.EvidenceItem;
import com.hartwig.hmftools.common.clinical.PatientTumorLocation;
import com.hartwig.hmftools.common.variant.structural.linx.LinxBreakend;
import com.hartwig.hmftools.common.variant.structural.linx.LinxFusion;
import com.hartwig.hmftools.patientreporter.actionability.ReportableEvidenceItemFactory;

import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public final class SvAnalyzer {

    private SvAnalyzer() {
    }

    @NotNull
    public static SvAnalysis run(@NotNull List<LinxFusion> fusions, @NotNull List<LinxBreakend> disruptions,
            @NotNull ActionabilityAnalyzer actionabilityAnalyzer, @Nullable PatientTumorLocation patientTumorLocation) {
        List<ReportableGeneDisruption> reportableGeneDisruptions = ReportableGeneDisruptionFactory.convert(disruptions);

        String primaryTumorLocation = patientTumorLocation != null ? patientTumorLocation.primaryTumorLocation() : null;

        Map<LinxFusion, List<EvidenceItem>> evidencePerFusion =
                actionabilityAnalyzer.evidenceForFusions(fusions, primaryTumorLocation);

        List<EvidenceItem> filteredEvidence = ReportableEvidenceItemFactory.toReportableFlatList(evidencePerFusion);

        return ImmutableSvAnalysis.builder()
                .reportableFusions(fusions)
                .reportableDisruptions(reportableGeneDisruptions)
                .evidenceItems(filteredEvidence)
                .build();
    }
}
