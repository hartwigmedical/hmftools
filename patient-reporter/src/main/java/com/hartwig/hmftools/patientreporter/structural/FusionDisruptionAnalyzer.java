package com.hartwig.hmftools.patientreporter.structural;

import java.io.IOException;
import java.util.List;
import java.util.Map;
import java.util.Objects;
import java.util.Set;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.actionability.ActionabilityAnalyzer;
import com.hartwig.hmftools.common.actionability.EvidenceItem;
import com.hartwig.hmftools.common.actionability.cancertype.CancerTypeMappingReading;
import com.hartwig.hmftools.common.ecrf.projections.PatientTumorLocation;
import com.hartwig.hmftools.common.purple.gene.GeneCopyNumber;
import com.hartwig.hmftools.common.variant.structural.annotation.GeneFusion;
import com.hartwig.hmftools.svannotation.analysis.StructuralVariantAnalysis;

import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public final class FusionDisruptionAnalyzer {

    private FusionDisruptionAnalyzer() {
    }

    @NotNull
    public static FusionDisruptionAnalysis run(@NotNull StructuralVariantAnalysis structuralVariantAnalysis,
            @NotNull List<GeneCopyNumber> exomeGeneCopyNumbers, @Nullable PatientTumorLocation patientTumorLocation,
            @NotNull ActionabilityAnalyzer actionabilityAnalyzer) throws IOException {
        List<ReportableGeneDisruption> reportableDisruptions = ReportableGeneDisruptionFactory.toReportableGeneDisruptions(
                structuralVariantAnalysis.reportableDisruptions(),
                exomeGeneCopyNumbers);

        List<GeneFusion> reportableFusions = structuralVariantAnalysis.reportableFusions();

        String primaryTumorLocation = patientTumorLocation != null ? patientTumorLocation.primaryTumorLocation() : Strings.EMPTY;
        CancerTypeMappingReading cancerTypeMappingReading = CancerTypeMappingReading.readingFile();
        String doidsPrimaryTumorLocation = cancerTypeMappingReading.doidsForPrimaryTumorLocation(primaryTumorLocation);

        Map<GeneFusion, List<EvidenceItem>> evidencePerFusion =
                findEvidenceForFusions(structuralVariantAnalysis.fusions(), doidsPrimaryTumorLocation, actionabilityAnalyzer);

        // KODU: Adding fusions with evidence if they are not selected yet.
        // KODU: Should reuse "favor canonical" rules from SV analyser here but have re-implemented for now.
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
            if (fusion.upstreamLinkedAnnotation().geneName().equals(newFusion.upstreamLinkedAnnotation().geneName()) &&
                    fusion.downstreamLinkedAnnotation().geneName().equals(newFusion.downstreamLinkedAnnotation().geneName())) {
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

    @NotNull
    private static Map<GeneFusion, List<EvidenceItem>> findEvidenceForFusions(@NotNull List<GeneFusion> fusions,
            @Nullable String doidsPrimaryTumorLocation, @NotNull ActionabilityAnalyzer actionabilityAnalyzer) {
        Set<String> actionableGenesInFusions = actionabilityAnalyzer.fusionAnalyzer().actionableGenes();

        List<GeneFusion> fusionsOnActionableGenes = fusions.stream()
                .filter(fusion -> actionableGenesInFusions.contains(fusion.upstreamLinkedAnnotation().geneName())
                        || actionableGenesInFusions.contains(fusion.downstreamLinkedAnnotation().geneName()))
                .collect(Collectors.toList());

        Map<GeneFusion, List<EvidenceItem>> evidenceItemsFusions = Maps.newHashMap();
        for (GeneFusion fusion : uniqueGeneFusions(fusionsOnActionableGenes)) {
            evidenceItemsFusions.put(fusion,
                    actionabilityAnalyzer.fusionAnalyzer()
                            .actionableFusions(doidsPrimaryTumorLocation, actionabilityAnalyzer.cancerTypeAnalyzer(), fusion));
        }

        return evidenceItemsFusions;
    }

    @NotNull
    private static List<GeneFusion> uniqueGeneFusions(@NotNull List<GeneFusion> fusions) {
        List<GeneFusion> allUniqueGeneFusions = Lists.newArrayList();
        Map<FiveThreePair, List<GeneFusion>> fusionsPerFiveThreePair = Maps.newHashMap();
        for (GeneFusion fusion : fusions) {
            FiveThreePair key = new FiveThreePair(fusion.upstreamLinkedAnnotation().geneName(),
                    fusion.downstreamLinkedAnnotation().geneName());
            List<GeneFusion> fusionsForKey = fusionsPerFiveThreePair.get(key);
            if (fusionsForKey == null) {
                fusionsForKey = Lists.newArrayList();
            }
            fusionsForKey.add(fusion);
            fusionsPerFiveThreePair.put(key, fusionsForKey);
        }

        for (Map.Entry<FiveThreePair, List<GeneFusion>> fusionsPerKey : fusionsPerFiveThreePair.entrySet()) {
            allUniqueGeneFusions.add(favorCanonical(fusionsPerKey.getValue()));
        }

        return allUniqueGeneFusions;
    }

    @NotNull
    private static GeneFusion favorCanonical(@NotNull List<GeneFusion> fusions) {
        for (GeneFusion fusion : fusions) {
            if (fusion.upstreamLinkedAnnotation().isCanonical() && fusion.downstreamLinkedAnnotation().isCanonical()) {
                return fusion;
            }
        }

        // KODU: If there is no canonical-canonical fusion, return the first one arbitrarily.
        assert !fusions.isEmpty();
        return fusions.get(0);
    }

    private static class FiveThreePair {

        @NotNull
        private final String fiveGene;
        @NotNull
        private final String threeGene;

        private FiveThreePair(@NotNull final String fiveGene, @NotNull final String threeGene) {
            this.fiveGene = fiveGene;
            this.threeGene = threeGene;
        }

        @Override
        public boolean equals(final Object o) {
            if (this == o) {
                return true;
            }
            if (o == null || getClass() != o.getClass()) {
                return false;
            }
            final FiveThreePair that = (FiveThreePair) o;
            return Objects.equals(fiveGene, that.fiveGene) && Objects.equals(threeGene, that.threeGene);
        }

        @Override
        public int hashCode() {
            return Objects.hash(fiveGene, threeGene);
        }
    }
}
