package com.hartwig.hmftools.patientreporter.structural;

import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.actionability.ActionabilityAnalyzer;
import com.hartwig.hmftools.common.actionability.EvidenceItem;
import com.hartwig.hmftools.common.ecrf.projections.PatientTumorLocation;
import com.hartwig.hmftools.common.purple.gene.GeneCopyNumber;
import com.hartwig.hmftools.common.variant.structural.annotation.GeneFusion;
import com.hartwig.hmftools.svannotation.analysis.StructuralVariantAnalysis;

import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public final class FusionDisruptionAnalyzer {

    private FusionDisruptionAnalyzer() {
    }

    @NotNull
    public static FusionDisruptionAnalysis run(@NotNull StructuralVariantAnalysis structuralVariantAnalysis,
            @NotNull List<GeneCopyNumber> exomeGeneCopyNumbers, @NotNull ActionabilityAnalyzer actionabilityAnalyzer,
            @Nullable PatientTumorLocation patientTumorLocation) {
        List<ReportableGeneDisruption> reportableDisruptions = ReportableGeneDisruptionFactory.toReportableGeneDisruptions(
                structuralVariantAnalysis.reportableDisruptions(),
                exomeGeneCopyNumbers);

        List<GeneFusion> reportableFusions = structuralVariantAnalysis.reportableFusions();

        String primaryTumorLocation = patientTumorLocation != null ? patientTumorLocation.primaryTumorLocation() : null;
        Map<GeneFusion, List<EvidenceItem>> evidencePerFusion =
                actionabilityAnalyzer.evidenceForFusions(structuralVariantAnalysis.fusions(), primaryTumorLocation);

        for (GeneFusion actionableFusion : evidencePerFusion.keySet()) {
            if (!evidencePerFusion.get(actionableFusion).isEmpty() && !fiveThreeCombinationExists(reportableFusions, actionableFusion)) {
                reportableFusions.add(actionableFusion);
            }
        }

        return ImmutableFusionDisruptionAnalysis.builder()
                .reportableFusions(ReportableGeneFusionFactory.toReportableGeneFusions(reportableFusions))
                .reportableDisruptions(reportableDisruptions)
                .evidenceItems(toList(evidencePerFusion))
                .build();
    }

    private static boolean fiveThreeCombinationExists(@NotNull List<GeneFusion> fusions, @NotNull GeneFusion newFusion) {
        for (GeneFusion fusion : fusions) {
            if (fusion.upstreamLinkedAnnotation().geneName().equals(newFusion.upstreamLinkedAnnotation().geneName())
                    && fusion.downstreamLinkedAnnotation().geneName().equals(newFusion.downstreamLinkedAnnotation().geneName())) {
                return true;
            }
        }

        return false;
    }

    @NotNull
    private static List<EvidenceItem> toList(@NotNull Map<GeneFusion, List<EvidenceItem>> evidencePerFusion) {
        List<EvidenceItem> evidenceItemList = Lists.newArrayList();
        for (List<EvidenceItem> items : evidencePerFusion.values()) {
            evidenceItemList.addAll(items);
        }
        return evidenceItemList;
    }
}
