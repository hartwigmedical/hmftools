package com.hartwig.hmftools.patientreporter.variants;

import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.stream.Collectors;

import com.google.common.collect.Maps;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.actionability.ActionabilityAnalyzer;
import com.hartwig.hmftools.common.actionability.EvidenceItem;
import com.hartwig.hmftools.common.purple.gene.GeneCopyNumber;
import com.hartwig.hmftools.common.variant.CodingEffect;
import com.hartwig.hmftools.common.variant.SomaticVariant;

import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public final class ActionabilityVariantAnalyzer {

    private static final Set<CodingEffect> CODING_EFFECTS =
            Sets.newHashSet(CodingEffect.SPLICE, CodingEffect.NONSENSE_OR_FRAMESHIFT, CodingEffect.MISSENSE);

    private ActionabilityVariantAnalyzer() {
    }

    @NotNull
    public static <T extends SomaticVariant> Map<T, List<EvidenceItem>> findEvidenceForVariants(
            @NotNull Set<String> actionableGenesVariants, @NotNull List<T> variants, @Nullable String doidsPrimaryTumorLocation,
            @NotNull ActionabilityAnalyzer actionabilityAnalyzerData) {
        Map<T, List<EvidenceItem>> evidenceItemsPerVariant = Maps.newHashMap();

        List<T> codingVariantsOnActionableGenesRanges = variants.stream()
                .filter(variant -> actionableGenesVariants.contains(variant.gene())
                        && CODING_EFFECTS.contains(variant.canonicalCodingEffect()))
                .collect(Collectors.toList());

        for (T variant : codingVariantsOnActionableGenesRanges) {
            evidenceItemsPerVariant.put(variant,
                    actionabilityAnalyzerData.variantAnalyzer()
                            .evidenceForSomaticVariant(variant, doidsPrimaryTumorLocation, actionabilityAnalyzerData.cancerTypeAnalyzer()));
        }
        return evidenceItemsPerVariant;
    }

    @NotNull
    public static Map<GeneCopyNumber, List<EvidenceItem>> detectCNVs(
            @NotNull Set<String> actionableGenesVariantsCNVs, @NotNull List<GeneCopyNumber> CNVs, @Nullable String doidsPrimaryTumorLocation,
            @NotNull ActionabilityAnalyzer actionabilityAnalyzerData) {
        Map<GeneCopyNumber, List<EvidenceItem>> evidenceItemsPerVariantCNVs = Maps.newHashMap();

        List<GeneCopyNumber> cnvsOnActionableGenes = CNVs.stream()
                .filter(geneCopyNumber -> actionableGenesVariantsCNVs.contains(geneCopyNumber.gene()))
                .collect(Collectors.toList());

        //        for (T CNV : variantsOnActionableGenes) {
        //            evidenceItemsPerVariantCNVs.put(CNV,
        //                    actionabilityAnalyzerData.cnvAnalyzer().evidenceForCopyNumberEvent(CNV, doidsPrimaryTumorLocation, actionabilityAnalyzerData));
        //        }

        return evidenceItemsPerVariantCNVs;
    }

    //    @NotNull
    //    public static <T extends GeneFusion> Map<T, FusionEvidenceItems> detectFusions(@NotNull List<T> fusion,
    //            @Nullable PatientTumorLocation patientTumorLocation) throws IOException {
    //        Map<T, FusionEvidenceItems> evidenceItemsPerVariantFusion = Maps.newHashMap();
    //        final String primaryTumorLocation = patientTumorLocation != null ? patientTumorLocation.primaryTumorLocation() : Strings.EMPTY;
    //        CancerTypeMappingReading cancerTypeMappingReading = CancerTypeMappingReading.readingFile();
    //        String doidsPrimaryTumorLocation = cancerTypeMappingReading.doidsForPrimaryTumorLocation(primaryTumorLocation);
    //        if (Files.exists(new File(FILE_ACTIONABILITY_FUSIONPAIRS).toPath())
    //                && Files.exists(new File(FILE_ACTIONABILITY_PROMISCUOUS_FIVE).toPath()) && Files.exists(new File(
    //                FILE_ACTIONABILITY_PROMISCUOUS_THREE).toPath()) && Files.exists(new File(FILE_CANCER_TUMORS_WITH_DOID).toPath())) {
    //            ActionabilityFusionAnalyzer analyzerFusions = ActionabilityFusionAnalyzer.loadFromFileFusions(FILE_ACTIONABILITY_FUSIONPAIRS,
    //                    FILE_ACTIONABILITY_PROMISCUOUS_FIVE,
    //                    FILE_ACTIONABILITY_PROMISCUOUS_THREE);
    //            CancerTypeAnalyzer cancerTypeAnalyzer = CancerTypeAnalyzer.loadFromFile(FILE_CANCER_TUMORS_WITH_DOID);
    //
    ////            Set<String> actionableGenesFusions = analyzerFusions.actionableGenes();
    ////            List<T> variantsOnActionableFusions = fusion.stream()
    ////                    .filter(fusionItem -> actionableGenesFusions.contains(fusionItem.reportable()))
    ////                    .collect(Collectors.toList());
    //
    ////            for (T fusionGene : variantsOnActionableFusions) {
    ////                evidenceItemsPerVariantFusion.put(fusionGene, analyzerFusions.actionableFusions(cancerTypeAnalyzer, doidsPrimaryTumorLocation));
    ////            }
    //
    //
    //        } else if (!Files.exists(new File(FILE_ACTIONABILITY_FUSIONPAIRS).toPath())) {
    //            LOGGER.warn("File does not exist: " + FILE_ACTIONABILITY_FUSIONPAIRS);
    //        } else if (!Files.exists(new File(FILE_ACTIONABILITY_PROMISCUOUS_FIVE).toPath())) {
    //            LOGGER.warn("File does not exist: " + FILE_ACTIONABILITY_PROMISCUOUS_FIVE);
    //        } else if (!Files.exists(new File(FILE_ACTIONABILITY_PROMISCUOUS_THREE).toPath())) {
    //            LOGGER.warn("File does not exist: " + FILE_ACTIONABILITY_PROMISCUOUS_THREE);
    //        } else if (!Files.exists(new File(FILE_CANCER_TUMORS_WITH_DOID).toPath())) {
    //            LOGGER.warn("File does not exist: " + FILE_CANCER_TUMORS_WITH_DOID);
    //        }
    //        return evidenceItemsPerVariantFusion;
    //    }
}
