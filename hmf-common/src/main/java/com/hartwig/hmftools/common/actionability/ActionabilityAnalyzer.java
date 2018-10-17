package com.hartwig.hmftools.common.actionability;

import java.io.File;
import java.io.IOException;

import com.hartwig.hmftools.common.actionability.cancertype.CancerTypeAnalyzer;
import com.hartwig.hmftools.common.actionability.cnv.ActionabilityCNVsAnalyzer;
import com.hartwig.hmftools.common.actionability.fusion.ActionabilityFusionAnalyzer;
import com.hartwig.hmftools.common.actionability.somaticvariant.SomaticVariantEvidenceAnalyzer;

import org.jetbrains.annotations.NotNull;

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
    private final ActionabilityCNVsAnalyzer cnvAnalyzer;
    @NotNull
    private final ActionabilityFusionAnalyzer fusionAnalyzer;
    @NotNull
    private final CancerTypeAnalyzer cancerTypeAnalyzer;

    @NotNull
    public static ActionabilityAnalyzer fromKnowledgebase(@NotNull String knowledgebasePath) throws IOException {
        String basePath = knowledgebasePath + File.separator;
        SomaticVariantEvidenceAnalyzer variantAnalyzer =
                SomaticVariantEvidenceAnalyzer.loadFromFileVariantsAndFileRanges(basePath + ACTIONABLE_VARIANT_FILE,
                        basePath + ACTIONABLE_RANGES_FILE);

        ActionabilityCNVsAnalyzer cnvAnalyzer = ActionabilityCNVsAnalyzer.loadFromFileCNVs(basePath + ACTIONABLE_CNV_FILE);

        ActionabilityFusionAnalyzer fusionAnalyzer = ActionabilityFusionAnalyzer.loadFromFileFusions(basePath + ACTIONABLE_FUSION_PAIR_FILE,
                basePath + ACTIONABLE_PROMISCUOUS_FIVE_FILE,
                basePath + ACTIONABLE_PROMISCUOUS_THREE_FILE);

        CancerTypeAnalyzer cancerTypeAnalyzer =
                CancerTypeAnalyzer.loadFromFile(knowledgebasePath + File.separator + CANCER_TYPE_DOID_MAPPING_FILE);

        return new ActionabilityAnalyzer(variantAnalyzer, cnvAnalyzer, fusionAnalyzer, cancerTypeAnalyzer);
    }

    private ActionabilityAnalyzer(@NotNull final SomaticVariantEvidenceAnalyzer variantAnalyzer,
            @NotNull final ActionabilityCNVsAnalyzer cnvAnalyzer,
            @NotNull final ActionabilityFusionAnalyzer fusionAnalyzer, @NotNull final CancerTypeAnalyzer cancerTypeAnalyzer) {
        this.variantAnalyzer = variantAnalyzer;
        this.cnvAnalyzer = cnvAnalyzer;
        this.fusionAnalyzer = fusionAnalyzer;
        this.cancerTypeAnalyzer = cancerTypeAnalyzer;
    }

    @NotNull
    public SomaticVariantEvidenceAnalyzer variantAnalyzer() {
        return variantAnalyzer;
    }

    @NotNull
    public ActionabilityCNVsAnalyzer cnvAnalyzer() {
        return cnvAnalyzer;
    }

    @NotNull
    public ActionabilityFusionAnalyzer fusionAnalyzer() {
        return fusionAnalyzer;
    }

    @NotNull
    public CancerTypeAnalyzer cancerTypeAnalyzer() {
        return cancerTypeAnalyzer;
    }
}
