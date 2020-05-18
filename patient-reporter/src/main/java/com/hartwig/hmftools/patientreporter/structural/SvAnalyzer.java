package com.hartwig.hmftools.patientreporter.structural;

import java.util.List;
import java.util.Map;

import com.hartwig.hmftools.common.actionability.ActionabilityAnalyzer;
import com.hartwig.hmftools.common.actionability.EvidenceItem;
import com.hartwig.hmftools.common.ecrf.projections.PatientTumorLocation;
import com.hartwig.hmftools.common.fusion.ReportableDisruption;
import com.hartwig.hmftools.common.fusion.ReportableGeneFusion;
import com.hartwig.hmftools.patientreporter.actionability.ReportableEvidenceItemFactory;

import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public final class SvAnalyzer {

    private SvAnalyzer() {
    }

    @NotNull
    public static SvAnalysis run(@NotNull List<ReportableGeneFusion> fusions, @NotNull List<ReportableDisruption> disruptions,
            @NotNull ActionabilityAnalyzer actionabilityAnalyzer, @Nullable PatientTumorLocation patientTumorLocation) {
        List<ReportableGeneDisruption> reportableGeneDisruptions = ReportableGeneDisruptionFactory.convert(disruptions);

        String primaryTumorLocation = patientTumorLocation != null ? patientTumorLocation.primaryTumorLocation() : null;

        Map<ReportableGeneFusion, List<EvidenceItem>> evidencePerFusion =
                actionabilityAnalyzer.evidenceForFusions(fusions, primaryTumorLocation);

        List<EvidenceItem> filteredEvidence = ReportableEvidenceItemFactory.toReportableFlatList(evidencePerFusion);

        return ImmutableSvAnalysis.builder()
                .reportableFusions(fusions)
                .reportableDisruptions(reportableGeneDisruptions)
                .evidenceItems(filteredEvidence)
                .build();
    }
}
