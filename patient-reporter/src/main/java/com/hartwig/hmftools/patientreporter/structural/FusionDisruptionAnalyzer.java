package com.hartwig.hmftools.patientreporter.structural;

import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.actionability.ActionabilityAnalyzer;
import com.hartwig.hmftools.common.actionability.EvidenceItem;
import com.hartwig.hmftools.common.ecrf.projections.PatientTumorLocation;
import com.hartwig.hmftools.common.variant.structural.annotation.ImmutableSimpleGeneFusion;
import com.hartwig.hmftools.common.variant.structural.annotation.SimpleGeneFusion;
import com.hartwig.hmftools.patientreporter.actionability.ReportableEvidenceItemFactory;

import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public final class FusionDisruptionAnalyzer {

    private FusionDisruptionAnalyzer() {
    }

    @NotNull
    public static FusionDisruptionAnalysis run(@NotNull List<ReportableGeneFusion> geneFusionsToReport,
            @NotNull List<ReportableGeneDisruption> geneDisruptionsToReport, @NotNull ActionabilityAnalyzer actionabilityAnalyzer,
            @Nullable PatientTumorLocation patientTumorLocation) {
        String primaryTumorLocation = patientTumorLocation != null ? patientTumorLocation.primaryTumorLocation() : null;
        Map<SimpleGeneFusion, List<EvidenceItem>> evidencePerFusion =
                actionabilityAnalyzer.evidenceForFusions(toSimpleGeneFusions(geneFusionsToReport), primaryTumorLocation);

        List<EvidenceItem> filteredEvidence = ReportableEvidenceItemFactory.reportableFlatList(evidencePerFusion);

        return ImmutableFusionDisruptionAnalysis.builder()
                .reportableFusions(geneFusionsToReport)
                .reportableDisruptions(geneDisruptionsToReport)
                .evidenceItems(filteredEvidence)
                .build();
    }

    @NotNull
    private static List<SimpleGeneFusion> toSimpleGeneFusions(@NotNull List<ReportableGeneFusion> fusions) {
        List<SimpleGeneFusion> simpleGeneFusions = Lists.newArrayList();
        for (ReportableGeneFusion fusionReport : fusions) {
            simpleGeneFusions.add(ImmutableSimpleGeneFusion.builder()
                    .fiveGene(fusionReport.geneStart())
                    .threeGene(fusionReport.geneEnd())
                    .build());
        }
        return simpleGeneFusions;
    }
}
