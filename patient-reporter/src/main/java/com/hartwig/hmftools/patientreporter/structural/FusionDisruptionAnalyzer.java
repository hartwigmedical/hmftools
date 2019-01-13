package com.hartwig.hmftools.patientreporter.structural;

import java.util.Collections;
import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.actionability.ActionabilityAnalyzer;
import com.hartwig.hmftools.common.actionability.EvidenceItem;
import com.hartwig.hmftools.common.ecrf.projections.PatientTumorLocation;
import com.hartwig.hmftools.common.purple.gene.GeneCopyNumber;
import com.hartwig.hmftools.common.variant.structural.annotation.GeneFusion;
import com.hartwig.hmftools.common.variant.structural.annotation.ImmutableSimpleGeneFusion;
import com.hartwig.hmftools.common.variant.structural.annotation.SimpleGeneFusion;
import com.hartwig.hmftools.common.variant.structural.annotation.StructuralVariantAnalysis;
import com.hartwig.hmftools.patientreporter.actionability.ReportableEvidenceItemFactory;

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
        Map<SimpleGeneFusion, List<EvidenceItem>> evidencePerFusion =
                actionabilityAnalyzer.evidenceForFusions(toSimpleGeneFusions(structuralVariantAnalysis.fusions()), primaryTumorLocation);

        List<EvidenceItem> filteredEvidence = ReportableEvidenceItemFactory.reportableFlatList(evidencePerFusion);

        // Add all fusions with filtered evidence that have not previously been added.
        for (Map.Entry<SimpleGeneFusion, List<EvidenceItem>> entry : evidencePerFusion.entrySet()) {
            SimpleGeneFusion fusion = entry.getKey();
            if (!fiveThreeCombinationExists(reportableFusions, fusion) && !Collections.disjoint(entry.getValue(), filteredEvidence)) {
                reportableFusions.add(pickFusionToReport(structuralVariantAnalysis.fusions(), fusion));
            }
        }

        return ImmutableFusionDisruptionAnalysis.builder()
                .reportableFusions(ReportableGeneFusionFactory.toReportableGeneFusions(reportableFusions))
                .reportableDisruptions(reportableDisruptions)
                .evidenceItems(filteredEvidence)
                .build();
    }

    private static boolean fiveThreeCombinationExists(@NotNull List<GeneFusion> fusions, @NotNull SimpleGeneFusion newFusion) {
        for (GeneFusion fusion : fusions) {
            if (fusion.upstreamTrans().geneName().equals(newFusion.fiveGene()) && fusion.downstreamTrans()
                    .geneName()
                    .equals(newFusion.threeGene())) {
                return true;
            }
        }

        return false;
    }

    @NotNull
    private static List<SimpleGeneFusion> toSimpleGeneFusions(@NotNull List<GeneFusion> fusions) {
        List<SimpleGeneFusion> simpleGeneFusions = Lists.newArrayList();
        for (GeneFusion fusionReport : fusions) {
            simpleGeneFusions.add(ImmutableSimpleGeneFusion.builder()
                    .fiveGene(fusionReport.upstreamTrans().geneName())
                    .threeGene(fusionReport.downstreamTrans().geneName())
                    .build());
        }
        return simpleGeneFusions;
    }

    @NotNull
    private static GeneFusion pickFusionToReport(@NotNull List<GeneFusion> fusions, @NotNull SimpleGeneFusion simpleFusion) {
        List<GeneFusion> fusionsMatchWithSimpleFusion = Lists.newArrayList();
        for (GeneFusion fusion : fusions) {
            if (fusion.upstreamTrans().geneName().equals(simpleFusion.fiveGene()) && fusion.downstreamTrans()
                    .geneName()
                    .equals(simpleFusion.threeGene())) {
                fusionsMatchWithSimpleFusion.add(fusion);
            }
        }

        for (GeneFusion fusion : fusionsMatchWithSimpleFusion) {
            if (fusion.upstreamTrans().isCanonical() && fusion.downstreamTrans().isCanonical()) {
                return fusion;
            }
        }

        // If there is no canonical-canonical fusion, return the first one arbitrarily.
        assert !fusions.isEmpty();
        return fusions.get(0);
    }
}
