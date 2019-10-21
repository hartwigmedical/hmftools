package com.hartwig.hmftools.common.actionability;

import java.io.File;
import java.io.IOException;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.stream.Collectors;

import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.actionability.cancertype.CancerTypeAnalyzer;
import com.hartwig.hmftools.common.actionability.cnv.CopyNumberEvidenceAnalyzer;
import com.hartwig.hmftools.common.actionability.cnv.CopyNumberEvidenceAnalyzerFactory;
import com.hartwig.hmftools.common.actionability.cnv.SignificantGeneCopyNumberFilter;
import com.hartwig.hmftools.common.actionability.fusion.FusionEvidenceAnalyzer;
import com.hartwig.hmftools.common.actionability.fusion.FusionEvidenceAnalyzerFactory;
import com.hartwig.hmftools.common.actionability.somaticvariant.SomaticVariantEvidenceAnalyzer;
import com.hartwig.hmftools.common.actionability.somaticvariant.SomaticVariantEvidenceAnalyzerFactory;
import com.hartwig.hmftools.common.purple.gene.GeneCopyNumber;
import com.hartwig.hmftools.common.variant.ReportableVariant;
import com.hartwig.hmftools.common.variant.SomaticVariant;
import com.hartwig.hmftools.common.variant.structural.annotation.ReportableGeneFusion;

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
    public static ActionabilityAnalyzer fromKnowledgebase(@NotNull String knowledgebaseDirectory) throws IOException {
        String basePath = knowledgebaseDirectory + File.separator;
        SomaticVariantEvidenceAnalyzer variantAnalyzer =
                SomaticVariantEvidenceAnalyzerFactory.loadFromFileVariantsAndFileRanges(basePath + ACTIONABLE_VARIANT_FILE,
                        basePath + ACTIONABLE_RANGES_FILE);

        CopyNumberEvidenceAnalyzer cnvAnalyzer = CopyNumberEvidenceAnalyzerFactory.loadFromFileCNVs(basePath + ACTIONABLE_CNV_FILE);

        FusionEvidenceAnalyzer fusionAnalyzer = FusionEvidenceAnalyzerFactory.loadFromFileFusions(basePath + ACTIONABLE_FUSION_PAIR_FILE,
                basePath + ACTIONABLE_PROMISCUOUS_FIVE_FILE,
                basePath + ACTIONABLE_PROMISCUOUS_THREE_FILE);

        CancerTypeAnalyzer cancerTypeAnalyzer = CancerTypeAnalyzer.createFromKnowledgeBase(basePath + CANCER_TYPE_DOID_MAPPING_FILE);

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

    @NotNull
    public SomaticVariantEvidenceAnalyzer variantAnalyzer() {
        return variantAnalyzer;
    }

    @NotNull
    public CopyNumberEvidenceAnalyzer cnvAnalyzer() {
        return copyNumberAnalyzer;
    }

    @NotNull
    public FusionEvidenceAnalyzer fusionAnalyzer() {
        return fusionAnalyzer;
    }

    @NotNull
    public Map<SomaticVariant, List<EvidenceItem>> evidenceForSomaticVariants(@NotNull List<SomaticVariant> variants,
            @Nullable String primaryTumorLocation) {
        Map<SomaticVariant, List<EvidenceItem>> evidencePerVariant = Maps.newHashMap();

        List<SomaticVariant> variantsOnActionableGenes = variants.stream()
                .filter(variant -> variantAnalyzer.actionableGenes().contains(variant.gene()))
                .collect(Collectors.toList());

        for (SomaticVariant variant : variantsOnActionableGenes) {
            evidencePerVariant.put(variant, variantAnalyzer.evidenceForSomaticVariant(variant, primaryTumorLocation, cancerTypeAnalyzer));
        }

        return evidencePerVariant;
    }

    @NotNull
    public Map<ReportableVariant, List<EvidenceItem>> evidenceForAllVariants(@NotNull List<ReportableVariant> variants,
            @Nullable String primaryTumorLocation) {
        Map<ReportableVariant, List<EvidenceItem>> evidencePerVariant = Maps.newHashMap();

        List<ReportableVariant> variantsOnActionableGenes = variants.stream()
                .filter(variant -> variantAnalyzer.actionableGenes().contains(variant.gene()))
                .collect(Collectors.toList());

        for (ReportableVariant variant : variantsOnActionableGenes) {
            evidencePerVariant.put(variant, variantAnalyzer.evidenceForAllVariant(variant, primaryTumorLocation, cancerTypeAnalyzer));
        }

        return evidencePerVariant;
    }

    @NotNull
    public Map<GeneCopyNumber, List<EvidenceItem>> evidenceForCopyNumbers(@NotNull List<GeneCopyNumber> geneCopyNumbers,
            @Nullable String primaryTumorLocation, double averageTumorPloidy) {
        Map<GeneCopyNumber, List<EvidenceItem>> evidencePerCopyNumber = Maps.newHashMap();

        List<GeneCopyNumber> geneCopyNumbersOnActionableGenes = geneCopyNumbers.stream()
                .filter(geneCopyNumber -> copyNumberAnalyzer.actionableGenes().contains(geneCopyNumber.gene()))
                .collect(Collectors.toList());

        // Filtering on significant events is not necessary but just to avoid unnecessary keys with empty evidence
        List<GeneCopyNumber> significantGeneCopyNumbersOnActionableGenes =
                SignificantGeneCopyNumberFilter.filterForSignificance(geneCopyNumbersOnActionableGenes, averageTumorPloidy);

        for (GeneCopyNumber geneCopyNumber : significantGeneCopyNumbersOnActionableGenes) {
            evidencePerCopyNumber.put(geneCopyNumber,
                    copyNumberAnalyzer.evidenceForCopyNumber(geneCopyNumber, averageTumorPloidy, primaryTumorLocation, cancerTypeAnalyzer));
        }

        return evidencePerCopyNumber;
    }

    @NotNull
    public Map<ReportableGeneFusion, List<EvidenceItem>> evidenceForFusions(@NotNull List<ReportableGeneFusion> fusions,
            @Nullable String primaryTumorLocation) {
        Map<ReportableGeneFusion, List<EvidenceItem>> evidencePerFusion = Maps.newHashMap();

        Set<ReportableGeneFusion> uniqueFusionsOnActionableGenes = fusions.stream()
                .filter(fusion -> fusionAnalyzer.actionableGenes().contains(fusion.geneStart()) || fusionAnalyzer.actionableGenes()
                        .contains(fusion.geneEnd()))
                .collect(Collectors.toSet());

        for (ReportableGeneFusion actionableFusion : uniqueFusionsOnActionableGenes) {
            evidencePerFusion.put(actionableFusion,
                    fusionAnalyzer.evidenceForFusion(actionableFusion, primaryTumorLocation, cancerTypeAnalyzer));
        }

        return evidencePerFusion;
    }
}
