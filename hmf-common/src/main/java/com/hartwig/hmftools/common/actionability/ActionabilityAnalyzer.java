package com.hartwig.hmftools.common.actionability;

import java.io.File;
import java.io.IOException;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.stream.Collectors;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.actionability.cancertype.CancerTypeAnalyzer;
import com.hartwig.hmftools.common.actionability.cnv.CopyNumberEvidenceAnalyzer;
import com.hartwig.hmftools.common.actionability.cnv.CopyNumberEvidenceAnalyzerFactory;
import com.hartwig.hmftools.common.actionability.fusion.FusionEvidenceAnalyzer;
import com.hartwig.hmftools.common.actionability.fusion.FusionEvidenceAnalyzerFactory;
import com.hartwig.hmftools.common.actionability.somaticvariant.SomaticVariantEvidenceAnalyzer;
import com.hartwig.hmftools.common.actionability.somaticvariant.SomaticVariantEvidenceAnalyzerFactory;
import com.hartwig.hmftools.common.purple.gene.GeneCopyNumber;
import com.hartwig.hmftools.common.variant.EnrichedSomaticVariant;
import com.hartwig.hmftools.common.variant.structural.annotation.SimpleGeneFusion;

import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public class ActionabilityAnalyzer {

    private static final String ACTIONABLE_VARIANT_FILE = "actionableVariants.tsv";
    private static final String ACTIONABLE_RANGES_FILE = "actionableRanges.tsv";
    private static final String ACTIONABLE_CNV_FILE = "actionableCNVs.tsv";
    private static final String ACTIONABLE_FUSION_PAIR_FILE = "actionableFusionPairs.tsv";
    private static final String ACTIONABLE_PROMISCUOUS_FIVE_FILE = "actionablePromiscuousFive.tsv";
    private static final String ACTIONABLE_PROMISCUOUS_THREE_FILE = "actionablePromiscuousThree.tsv";

    private static final String CANCER_TYPE_DOID_MAPPING_FILE = "knowledgebaseCancerTypes.tsv";

    @NotNull
    private final SomaticVariantEvidenceAnalyzer variantAnalyzer;
    @NotNull
    private final CopyNumberEvidenceAnalyzer copyNumberAnalyzer;
    @NotNull
    private final FusionEvidenceAnalyzer fusionAnalyzer;
    @NotNull
    private final CancerTypeAnalyzer cancerTypeAnalyzer;

    @NotNull
    public static ActionabilityAnalyzer fromKnowledgebase(@NotNull String knowledgebasePath) throws IOException {
        String basePath = knowledgebasePath + File.separator;
        SomaticVariantEvidenceAnalyzer variantAnalyzer =
                SomaticVariantEvidenceAnalyzerFactory.loadFromFileVariantsAndFileRanges(basePath + ACTIONABLE_VARIANT_FILE,
                        basePath + ACTIONABLE_RANGES_FILE);

        CopyNumberEvidenceAnalyzer cnvAnalyzer = CopyNumberEvidenceAnalyzerFactory.loadFromFileCNVs(basePath + ACTIONABLE_CNV_FILE);

        FusionEvidenceAnalyzer fusionAnalyzer = FusionEvidenceAnalyzerFactory.loadFromFileFusions(basePath + ACTIONABLE_FUSION_PAIR_FILE,
                basePath + ACTIONABLE_PROMISCUOUS_FIVE_FILE,
                basePath + ACTIONABLE_PROMISCUOUS_THREE_FILE);

        CancerTypeAnalyzer cancerTypeAnalyzer =
                CancerTypeAnalyzer.createFromKnowledgeBase(knowledgebasePath + File.separator + CANCER_TYPE_DOID_MAPPING_FILE);

        return new ActionabilityAnalyzer(variantAnalyzer, cnvAnalyzer, fusionAnalyzer, cancerTypeAnalyzer);
    }

    private ActionabilityAnalyzer(@NotNull final SomaticVariantEvidenceAnalyzer variantAnalyzer,
            @NotNull final CopyNumberEvidenceAnalyzer copyNumberAnalyzer, @NotNull final FusionEvidenceAnalyzer fusionAnalyzer,
            @NotNull final CancerTypeAnalyzer cancerTypeAnalyzer) {
        this.variantAnalyzer = variantAnalyzer;
        this.copyNumberAnalyzer = copyNumberAnalyzer;
        this.fusionAnalyzer = fusionAnalyzer;
        this.cancerTypeAnalyzer = cancerTypeAnalyzer;
    }

    @VisibleForTesting
    @NotNull
    SomaticVariantEvidenceAnalyzer variantAnalyzer() {
        return variantAnalyzer;
    }

    @VisibleForTesting
    @NotNull
    CopyNumberEvidenceAnalyzer cnvAnalyzer() {
        return copyNumberAnalyzer;
    }

    @VisibleForTesting
    @NotNull
    FusionEvidenceAnalyzer fusionAnalyzer() {
        return fusionAnalyzer;
    }

    @NotNull
    public Map<EnrichedSomaticVariant, List<EvidenceItem>> evidenceForSomaticVariants(@NotNull List<EnrichedSomaticVariant> variants,
            @Nullable String primaryTumorLocation) {
        Map<EnrichedSomaticVariant, List<EvidenceItem>> evidencePerVariant = Maps.newHashMap();

        List<EnrichedSomaticVariant> variantsOnActionableGenes = variants.stream()
                .filter(variant -> variantAnalyzer.actionableGenes().contains(variant.gene()))
                .collect(Collectors.toList());

        for (EnrichedSomaticVariant variant : variantsOnActionableGenes) {
            evidencePerVariant.put(variant, variantAnalyzer.evidenceForSomaticVariant(variant, primaryTumorLocation, cancerTypeAnalyzer));
        }

        return evidencePerVariant;
    }

    @NotNull
    public Map<GeneCopyNumber, List<EvidenceItem>> evidenceForCopyNumbers(@NotNull List<GeneCopyNumber> geneCopyNumbers,
            @Nullable String primaryTumorLocation) {
        Map<GeneCopyNumber, List<EvidenceItem>> evidencePerCopyNumber = Maps.newHashMap();

        // TODO (KODU): Should also filter on significant event rather than assume caller has filtered already.
        List<GeneCopyNumber> geneCopyNumbersOnActionableGenes = geneCopyNumbers.stream()
                .filter(geneCopyNumber -> copyNumberAnalyzer.actionableGenes().contains(geneCopyNumber.gene()))
                .collect(Collectors.toList());

        for (GeneCopyNumber geneCopyNumber : geneCopyNumbersOnActionableGenes) {
            evidencePerCopyNumber.put(geneCopyNumber,
                    copyNumberAnalyzer.evidenceForCopyNumber(geneCopyNumber, primaryTumorLocation, cancerTypeAnalyzer));
        }

        return evidencePerCopyNumber;
    }

    @NotNull
    public Map<SimpleGeneFusion, List<EvidenceItem>> evidenceForFusions(@NotNull List<SimpleGeneFusion> fusions,
            @Nullable String primaryTumorLocation) {
        Map<SimpleGeneFusion, List<EvidenceItem>> evidencePerFusion = Maps.newHashMap();

        Set<SimpleGeneFusion> uniqueFusionsOnActionableGenes = fusions.stream()
                .filter(fusion -> fusionAnalyzer.actionableGenes().contains(fusion.fiveGene()) || fusionAnalyzer.actionableGenes()
                        .contains(fusion.threeGene()))
                .collect(Collectors.toSet());

        for (SimpleGeneFusion actionableFusion : uniqueFusionsOnActionableGenes) {
            evidencePerFusion.put(actionableFusion,
                    fusionAnalyzer.evidenceForFusion(actionableFusion, primaryTumorLocation, cancerTypeAnalyzer));
        }

        return evidencePerFusion;
    }
}
