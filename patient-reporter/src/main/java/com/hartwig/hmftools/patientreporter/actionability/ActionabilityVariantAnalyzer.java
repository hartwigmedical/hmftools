package com.hartwig.hmftools.patientreporter.actionability;

import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.actionability.ActionabilityAnalyzer;
import com.hartwig.hmftools.common.actionability.EvidenceItem;
import com.hartwig.hmftools.common.purple.gene.GeneCopyNumber;
import com.hartwig.hmftools.common.variant.structural.annotation.GeneFusion;

import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public final class ActionabilityVariantAnalyzer {

    private ActionabilityVariantAnalyzer() {
    }

    @NotNull
    public static Map<GeneCopyNumber, List<EvidenceItem>> findEvidenceForCopyNumber(@NotNull Set<String> actionableGenesVariantsCNVs,
            @NotNull List<GeneCopyNumber> CNVs, @Nullable String doidsPrimaryTumorLocation,
            @NotNull ActionabilityAnalyzer actionabilityAnalyzerData, final double purplePloidy) {
        Map<GeneCopyNumber, List<EvidenceItem>> evidenceItemsCopyNumber = Maps.newHashMap();

        List<GeneCopyNumber> cnvsOnActionableGenes = CNVs.stream()
                .filter(geneCopyNumber -> actionableGenesVariantsCNVs.contains(geneCopyNumber.gene()))
                .collect(Collectors.toList());

        for (GeneCopyNumber CNV : cnvsOnActionableGenes) {
            evidenceItemsCopyNumber.put(CNV,
                    actionabilityAnalyzerData.cnvAnalyzer()
                            .evidenceForCopyNumberEvent(CNV,
                                    doidsPrimaryTumorLocation,
                                    actionabilityAnalyzerData.cancerTypeAnalyzer(),
                                    purplePloidy));
        }

        return evidenceItemsCopyNumber;
    }

    @NotNull
    public static Map<GeneFusion, List<EvidenceItem>> findEvidenceForFusions(@NotNull Set<String> actionableGenesFusions,
            @NotNull List<GeneFusion> fusions, @Nullable String doidsPrimaryTumorLocation,
            @NotNull ActionabilityAnalyzer actionabilityAnalyzerData) {
        Map<GeneFusion, List<EvidenceItem>> evidenceItemsFusions = Maps.newHashMap();
        List<GeneFusion> fusionsOnActionableGenes = fusions.stream()
                .filter(fusion -> actionableGenesFusions.contains(fusion.downstreamLinkedAnnotation().geneName())
                        || actionableGenesFusions.contains(fusion.upstreamLinkedAnnotation().geneName()))
                .collect(Collectors.toList());

        for (GeneFusion fusion : uniqueGeneFusions(fusionsOnActionableGenes)) {
            evidenceItemsFusions.put(fusion,
                    actionabilityAnalyzerData.fusionAnalyzer()
                            .actionableFusions(doidsPrimaryTumorLocation, actionabilityAnalyzerData.cancerTypeAnalyzer(), fusion));
        }

        return evidenceItemsFusions;
    }

    @NotNull
    private static List<GeneFusion> uniqueGeneFusions(@NotNull List<GeneFusion> fusions) {
        List<GeneFusion> allUniqueGeneFusions = Lists.newArrayList();
        List<String> uniqueGeneFive = Lists.newArrayList();
        List<String> uniqueGeneThree = Lists.newArrayList();
        for (GeneFusion fusion : fusions) {
            if (!uniqueGeneFive.contains(fusion.upstreamLinkedAnnotation().geneName())
                    || !uniqueGeneThree.contains(fusion.downstreamLinkedAnnotation().geneName())) {
                uniqueGeneFive.add(fusion.upstreamLinkedAnnotation().geneName());
                uniqueGeneThree.add(fusion.downstreamLinkedAnnotation().geneName());
                allUniqueGeneFusions.add(fusion);
            }
        }
        return allUniqueGeneFusions;
    }
}
