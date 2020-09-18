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
import com.hartwig.hmftools.common.actionability.variant.VariantEvidenceAnalyzer;
import com.hartwig.hmftools.common.actionability.variant.VariantEvidenceAnalyzerFactory;
import com.hartwig.hmftools.common.purple.copynumber.ReportableGainLoss;
import com.hartwig.hmftools.common.variant.Variant;
import com.hartwig.hmftools.common.fusion.ReportableGeneFusion;

import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public class ActionabilityAnalyzer {

    private static final String ACTIONABLE_VARIANT_TSV = "actionableVariants.tsv";
    private static final String ACTIONABLE_RANGES_TSV = "actionableRanges.tsv";
    private static final String ACTIONABLE_CNV_TSV = "actionableCNVs.tsv";
    private static final String ACTIONABLE_FUSION_PAIR_TSV = "actionableFusionPairs.tsv";
    private static final String ACTIONABLE_PROMISCUOUS_FIVE_TSV = "actionablePromiscuousFive.tsv";
    private static final String ACTIONABLE_PROMISCUOUS_THREE_TSV = "actionablePromiscuousThree.tsv";

    private static final String CANCER_TYPE_DOID_MAPPING_TSV = "knowledgebaseCancerTypes.tsv";

    @NotNull
    private final VariantEvidenceAnalyzer variantAnalyzer;
    @NotNull
    private final CopyNumberEvidenceAnalyzer copyNumberAnalyzer;
    @NotNull
    private final FusionEvidenceAnalyzer fusionAnalyzer;
    @NotNull
    private final CancerTypeAnalyzer cancerTypeAnalyzer;

    @NotNull
    public static ActionabilityAnalyzer fromKnowledgebase(@NotNull String knowledgebaseDirectory) throws IOException {
        String basePath = knowledgebaseDirectory + File.separator;
        VariantEvidenceAnalyzer variantAnalyzer =
                VariantEvidenceAnalyzerFactory.loadFromFileVariantsAndFileRanges(basePath + ACTIONABLE_VARIANT_TSV,
                        basePath + ACTIONABLE_RANGES_TSV);

        CopyNumberEvidenceAnalyzer cnvAnalyzer = CopyNumberEvidenceAnalyzerFactory.loadFromFileCNVs(basePath + ACTIONABLE_CNV_TSV);

        FusionEvidenceAnalyzer fusionAnalyzer = FusionEvidenceAnalyzerFactory.loadFromFileFusions(basePath + ACTIONABLE_FUSION_PAIR_TSV,
                basePath + ACTIONABLE_PROMISCUOUS_FIVE_TSV,
                basePath + ACTIONABLE_PROMISCUOUS_THREE_TSV);

        CancerTypeAnalyzer cancerTypeAnalyzer = CancerTypeAnalyzer.createFromKnowledgeBase(basePath + CANCER_TYPE_DOID_MAPPING_TSV);

        return new ActionabilityAnalyzer(variantAnalyzer, cnvAnalyzer, fusionAnalyzer, cancerTypeAnalyzer);
    }

    private ActionabilityAnalyzer(@NotNull final VariantEvidenceAnalyzer variantAnalyzer,
            @NotNull final CopyNumberEvidenceAnalyzer copyNumberAnalyzer, @NotNull final FusionEvidenceAnalyzer fusionAnalyzer,
            @NotNull final CancerTypeAnalyzer cancerTypeAnalyzer) {
        this.variantAnalyzer = variantAnalyzer;
        this.copyNumberAnalyzer = copyNumberAnalyzer;
        this.fusionAnalyzer = fusionAnalyzer;
        this.cancerTypeAnalyzer = cancerTypeAnalyzer;
    }

    @NotNull
    public VariantEvidenceAnalyzer variantAnalyzer() {
        return variantAnalyzer;
    }

    @NotNull
    public CopyNumberEvidenceAnalyzer cnvAnalyzer() {
        return copyNumberAnalyzer;
    }

    @NotNull
    @VisibleForTesting
    FusionEvidenceAnalyzer fusionAnalyzer() {
        return fusionAnalyzer;
    }

    @NotNull
    public Map<Variant, List<EvidenceItem>> evidenceForAllVariants(@NotNull List<? extends Variant> variants,
            @Nullable String primaryTumorLocation) {
        Map<Variant, List<EvidenceItem>> evidencePerVariant = Maps.newHashMap();

        List<Variant> variantsOnActionableGenes = variants.stream()
                .filter(variant -> variantAnalyzer.actionableGenes().contains(variant.gene()))
                .collect(Collectors.toList());

        for (Variant variant : variantsOnActionableGenes) {
            evidencePerVariant.put(variant, variantAnalyzer.evidenceForVariant(variant, primaryTumorLocation, cancerTypeAnalyzer));
        }

        return evidencePerVariant;
    }

    @NotNull
    public Map<ReportableGainLoss, List<EvidenceItem>> evidenceForCopyNumbers(@NotNull List<ReportableGainLoss> reportableGainsAndLosses,
            @Nullable String primaryTumorLocation, double averageTumorPloidy) {
        Map<ReportableGainLoss, List<EvidenceItem>> evidencePerCopyNumber = Maps.newHashMap();

        List<ReportableGainLoss> geneCopyNumbersOnActionableGenes = reportableGainsAndLosses.stream()
                .filter(reportableGainAndLoss -> copyNumberAnalyzer.actionableGenes().contains(reportableGainAndLoss.gene()))
                .collect(Collectors.toList());

        // Filtering on significant events is not necessary but just to avoid unnecessary keys with empty evidence
//        List<GeneCopyNumber> significantGeneCopyNumbersOnActionableGenes =
//                SignificantGeneCopyNumberFilter.filterForSignificance(geneCopyNumbersOnActionableGenes, averageTumorPloidy);

        for (ReportableGainLoss reportableGainLoss : geneCopyNumbersOnActionableGenes) {
            evidencePerCopyNumber.put(reportableGainLoss,
                    copyNumberAnalyzer.evidenceForCopyNumber(reportableGainLoss, averageTumorPloidy, primaryTumorLocation, cancerTypeAnalyzer));
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
