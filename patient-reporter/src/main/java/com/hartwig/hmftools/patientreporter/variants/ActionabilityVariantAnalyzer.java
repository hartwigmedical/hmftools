package com.hartwig.hmftools.patientreporter.variants;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.stream.Collectors;

import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.actionability.ActionabilityAnalyzer;
import com.hartwig.hmftools.common.actionability.cnv.ActionabilityCNVsEvidenceItems;
import com.hartwig.hmftools.common.actionability.somaticvariant.ActionabilityRangeEvidenceItem;
import com.hartwig.hmftools.common.actionability.somaticvariant.VariantEvidenceItems;
import com.hartwig.hmftools.common.purple.gene.GeneCopyNumber;
import com.hartwig.hmftools.common.variant.SomaticVariant;

import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public final class ActionabilityVariantAnalyzer {

    private ActionabilityVariantAnalyzer() {
    }

    @NotNull
    public static <T extends SomaticVariant> Map<T, VariantEvidenceItems> detectVariants(@NotNull Set<String> actionableGenesVariants,
            @NotNull List<T> variants, @Nullable String doidsPrimaryTumorLocation,
            @NotNull ActionabilityAnalyzer actionabilityAnalyzerData) {
        Map<T, VariantEvidenceItems> evidenceItemsPerVariant = Maps.newHashMap();

        List<T> variantsOnActionableGenes =
                variants.stream().filter(variant -> actionableGenesVariants.contains(variant.gene())).collect(Collectors.toList());

        for (T variant : variantsOnActionableGenes) {
            evidenceItemsPerVariant.put(variant,
                    actionabilityAnalyzerData.variantAnalyzer()
                            .actionableVariants(variant, doidsPrimaryTumorLocation, actionabilityAnalyzerData));
        }
        return evidenceItemsPerVariant;
    }

    @NotNull
    public static <T extends SomaticVariant> Map<T, ActionabilityRangeEvidenceItem> detectVariantsRanges(
            @NotNull Set<String> actionableGenesVariantsRanges, @NotNull List<T> variants, @Nullable String doidsPrimaryTumorLocation,
            @NotNull ActionabilityAnalyzer actionabilityAnalyzerData) {
        Map<T, ActionabilityRangeEvidenceItem> evidenceItemsPerVariantRange = Maps.newHashMap();

        List<T> variantsOnActionableGenesRanges =
                variants.stream().filter(variant -> actionableGenesVariantsRanges.contains(variant.gene())).collect(Collectors.toList());

        for (T variant : variantsOnActionableGenesRanges) {
            evidenceItemsPerVariantRange.put(variant,
                    actionabilityAnalyzerData.variantAnalyzer()
                            .actionableRange(variant, doidsPrimaryTumorLocation, actionabilityAnalyzerData));
        }

        return evidenceItemsPerVariantRange;
    }

    @NotNull
    public static <T extends GeneCopyNumber> Map<T, ActionabilityCNVsEvidenceItems> detectCNVs(
            @NotNull Set<String> actionableGenesVariantsCNVs, @NotNull List<T> CNVs, @Nullable String doidsPrimaryTumorLocation,
            @NotNull ActionabilityAnalyzer actionabilityAnalyzerData) {
        Map<T, ActionabilityCNVsEvidenceItems> evidenceItemsPerVariantCNVs = Maps.newHashMap();

        List<T> variantsOnActionableGenes =
                CNVs.stream().filter(geneCopyNumber -> actionableGenesVariantsCNVs.contains(geneCopyNumber.gene())).collect(Collectors.toList());

        for (T CNV : variantsOnActionableGenes) {
            evidenceItemsPerVariantCNVs.put(CNV, actionabilityAnalyzerData.cnvAnalyzer().actionableCNVs(CNV, doidsPrimaryTumorLocation, actionabilityAnalyzerData));
        }

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
