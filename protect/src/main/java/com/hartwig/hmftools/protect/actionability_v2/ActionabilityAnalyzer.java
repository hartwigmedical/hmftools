package com.hartwig.hmftools.protect.actionability_v2;

import java.io.File;
import java.io.IOException;

import com.hartwig.hmftools.protect.actionability.cancertype.CancerTypeAnalyzer;
import com.hartwig.hmftools.protect.actionability.cnv.CopyNumberEvidenceAnalyzer;
import com.hartwig.hmftools.protect.actionability_v2.fusion.FusionEvidenceAnalyzer;
import com.hartwig.hmftools.protect.actionability_v2.fusion.FusionEvidenceAnalyzerFactory;
import com.hartwig.hmftools.protect.actionability_v2.gene.GeneEvidenceAnalyzer;
import com.hartwig.hmftools.protect.actionability_v2.gene.GeneEvidenceAnalyzerFactory;
import com.hartwig.hmftools.protect.actionability_v2.range.RangeEvidenceAnalyzer;
import com.hartwig.hmftools.protect.actionability_v2.range.RangeEvidenceAnalyzerFactory;
import com.hartwig.hmftools.protect.actionability_v2.signature.SignatureEvidenceAnalyzer;
import com.hartwig.hmftools.protect.actionability_v2.signature.SignatureEvidenceAnalyzerFactory;
import com.hartwig.hmftools.protect.actionability_v2.variant.VariantEvidenceAnalyzer;
import com.hartwig.hmftools.protect.actionability_v2.variant.VariantEvidenceAnalyzerFactory;

import org.jetbrains.annotations.NotNull;

public class ActionabilityAnalyzer {

    private static final String ACTIONABLE_FUSION_TSV = "actionableFusion.tsv";
    private static final String ACTIONABLE_GENE_TSV = "actionableGene.tsv";
    private static final String ACTIONABLE_RANGE_TSV = "actionableRange.tsv";
    private static final String ACTIONABLE_SIGNATURE_TSV = "actionableSignature.tsv";
    private static final String ACTIONABLE_VARIANT_TSV = "actionableVariant.tsv";

    @NotNull
    private final FusionEvidenceAnalyzer fusionEvidenceAnalyzer;
    @NotNull
    private final GeneEvidenceAnalyzer geneEvidenceAnalyzer;
    @NotNull
    private final RangeEvidenceAnalyzer rangeEvidenceAnalyzer;
    @NotNull
    private final SignatureEvidenceAnalyzer signatureEvidenceAnalyzer;
    @NotNull
    private final VariantEvidenceAnalyzer variantEvidenceAnalyzer;

    @NotNull
    public static ActionabilityAnalyzer fromKnowledgebase(@NotNull String knowledgebaseDirectory) throws IOException {

        String basePath = knowledgebaseDirectory + File.separator;
        FusionEvidenceAnalyzer actionableFusion = FusionEvidenceAnalyzerFactory.loadFromFileFusion(basePath + ACTIONABLE_FUSION_TSV);
        GeneEvidenceAnalyzer actionableGene = GeneEvidenceAnalyzerFactory.loadFromFileGene(basePath + ACTIONABLE_GENE_TSV);
        RangeEvidenceAnalyzer actionableRange = RangeEvidenceAnalyzerFactory.loadFromFileRange(basePath + ACTIONABLE_RANGE_TSV);
        SignatureEvidenceAnalyzer actionableSignature =
                SignatureEvidenceAnalyzerFactory.loadFromFileSignature(basePath + ACTIONABLE_SIGNATURE_TSV);
        VariantEvidenceAnalyzer actionableVariant = VariantEvidenceAnalyzerFactory.loadFromFileVariant(basePath + ACTIONABLE_VARIANT_TSV);

        return new ActionabilityAnalyzer(actionableFusion, actionableGene, actionableRange, actionableSignature, actionableVariant);
    }

    private ActionabilityAnalyzer(@NotNull final FusionEvidenceAnalyzer actionableFusion,
            @NotNull final GeneEvidenceAnalyzer actionableGene, @NotNull final RangeEvidenceAnalyzer actionableRange,
            @NotNull final SignatureEvidenceAnalyzer actionableSignature, @NotNull final VariantEvidenceAnalyzer actionableVariant) {
        this.fusionEvidenceAnalyzer = actionableFusion;
        this.geneEvidenceAnalyzer = actionableGene;
        this.rangeEvidenceAnalyzer = actionableRange;
        this.signatureEvidenceAnalyzer = actionableSignature;
        this.variantEvidenceAnalyzer = actionableVariant;
    }
}
