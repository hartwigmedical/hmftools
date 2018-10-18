package com.hartwig.hmftools.patientreporter.variants;

import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.stream.Collectors;

import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.actionability.ActionabilityAnalyzer;
import com.hartwig.hmftools.common.actionability.EvidenceItem;
import com.hartwig.hmftools.common.purple.gene.GeneCopyNumber;
import com.hartwig.hmftools.common.variant.SomaticVariant;
import com.hartwig.hmftools.common.variant.structural.annotation.GeneFusion;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public final class ActionabilityVariantAnalyzer {
    private static final Logger LOGGER = LogManager.getLogger(ActionabilityVariantAnalyzer.class);

    private ActionabilityVariantAnalyzer() {
    }

    @NotNull
    public static <T extends SomaticVariant> Map<T, List<EvidenceItem>> findEvidenceForVariants(
            @NotNull Set<String> actionableGenesVariants, @NotNull List<T> variants, @Nullable String doidsPrimaryTumorLocation,
            @NotNull ActionabilityAnalyzer actionabilityAnalyzerData) {
        Map<T, List<EvidenceItem>> evidenceItemsPerVariant = Maps.newHashMap();
        List<T> variantsOnActionableGenesRanges =
                variants.stream().filter(variant -> actionableGenesVariants.contains(variant.gene())).collect(Collectors.toList());

        for (T variant : variantsOnActionableGenesRanges) {
            evidenceItemsPerVariant.put(variant,
                    actionabilityAnalyzerData.variantAnalyzer()
                            .evidenceForSomaticVariant(variant, doidsPrimaryTumorLocation, actionabilityAnalyzerData.cancerTypeAnalyzer()));
        }
        return evidenceItemsPerVariant;
    }

    @NotNull
    public static <T extends GeneCopyNumber> Map<T, List<EvidenceItem>> findEvidenceForCopyNumber(
            @NotNull Set<String> actionableGenesVariantsCNVs, @NotNull List<T> CNVs, @Nullable String doidsPrimaryTumorLocation,
            @NotNull ActionabilityAnalyzer actionabilityAnalyzerData, final double purplePloidy) {
        Map<T, List<EvidenceItem>> evidenceItemsCopyNumber = Maps.newHashMap();

        List<T> cnvsOnActionableGenes = CNVs.stream()
                .filter(geneCopyNumber -> actionableGenesVariantsCNVs.contains(geneCopyNumber.gene()))
                .collect(Collectors.toList());

        for (T CNV : cnvsOnActionableGenes) {
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
    public static <T extends GeneFusion> Map<T, List<EvidenceItem>> findEvidenceFusion(@NotNull Set<String> actionableGenesFusions,
            @NotNull List<T> fusions, @Nullable String doidsPrimaryTumorLocation,
            @NotNull ActionabilityAnalyzer actionabilityAnalyzerData) {
        Map<T, List<EvidenceItem>> evidenceItemsFusions = Maps.newHashMap();
        List<T> fusionsOnActionableGenes = fusions.stream()
                .filter(fusion -> actionableGenesFusions.contains(fusion.downstreamLinkedAnnotation().geneName())
                        && actionableGenesFusions.contains(fusion.upstreamLinkedAnnotation().geneName()))
                .collect(Collectors.toList());

        for (T fusion : fusionsOnActionableGenes) {
            evidenceItemsFusions.put(fusion,
                    actionabilityAnalyzerData.fusionAnalyzer()
                            .actionableFusions(doidsPrimaryTumorLocation, actionabilityAnalyzerData.cancerTypeAnalyzer(), fusion));
        }

        return evidenceItemsFusions;

    }
}
