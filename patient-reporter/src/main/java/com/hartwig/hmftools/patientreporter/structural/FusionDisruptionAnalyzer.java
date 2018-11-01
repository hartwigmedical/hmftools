package com.hartwig.hmftools.patientreporter.structural;

import java.util.Collections;
import java.util.List;
import java.util.Map;

import com.hartwig.hmftools.common.actionability.ActionabilityAnalyzer;
import com.hartwig.hmftools.common.actionability.EvidenceItem;
import com.hartwig.hmftools.common.ecrf.projections.PatientTumorLocation;
import com.hartwig.hmftools.common.purple.gene.GeneCopyNumber;
import com.hartwig.hmftools.common.variant.structural.annotation.GeneFusion;
import com.hartwig.hmftools.patientreporter.actionability.ReportableEvidenceItemFactory;
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

        List<EvidenceItem> filteredEvidence = ReportableEvidenceItemFactory.reportableFlatList(evidencePerFusion);

        // KODU: Add all fusions with filtered evidence that have not previously been added.
        for (Map.Entry<GeneFusion, List<EvidenceItem>> entry : evidencePerFusion.entrySet()) {
            GeneFusion fusion = entry.getKey();
            if (!fiveThreeCombinationExists(reportableFusions, fusion) && !Collections.disjoint(entry.getValue(), filteredEvidence)) {
                reportableFusions.add(fusion);
            }
        }

        return ImmutableFusionDisruptionAnalysis.builder()
                .reportableFusions(ReportableGeneFusionFactory.toReportableGeneFusions(reportableFusions))
                .reportableDisruptions(reportableDisruptions)
                .evidenceItems(filteredEvidence)
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
}
