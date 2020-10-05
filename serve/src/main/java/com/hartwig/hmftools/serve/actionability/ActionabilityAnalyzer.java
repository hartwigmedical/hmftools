package com.hartwig.hmftools.serve.actionability;

import java.io.File;
import java.io.IOException;

import com.google.common.annotations.VisibleForTesting;
import com.hartwig.hmftools.serve.actionability.fusion.FusionEvidenceAnalyzer;
import com.hartwig.hmftools.serve.actionability.fusion.FusionEvidenceAnalyzerFactory;
import com.hartwig.hmftools.serve.actionability.gene.GeneEvidenceAnalyzer;
import com.hartwig.hmftools.serve.actionability.gene.GeneEvidenceAnalyzerFactory;
import com.hartwig.hmftools.serve.actionability.range.RangeEvidenceAnalyzer;
import com.hartwig.hmftools.serve.actionability.range.RangeEvidenceAnalyzerFactory;
import com.hartwig.hmftools.serve.actionability.signature.SignatureEvidenceAnalyzer;
import com.hartwig.hmftools.serve.actionability.signature.SignatureEvidenceAnalyzerFactory;
import com.hartwig.hmftools.serve.actionability.variant.VariantEvidenceAnalyzer;
import com.hartwig.hmftools.serve.actionability.variant.VariantEvidenceAnalyzerFactory;

import org.jetbrains.annotations.NotNull;

public class ActionabilityAnalyzer {

    private static final String ACTIONABLE_FUSION_TSV = "actionableFusions.tsv";
    private static final String ACTIONABLE_GENE_TSV = "actionableGenes.tsv";
    private static final String ACTIONABLE_RANGE_TSV = "actionableRanges.tsv";
    private static final String ACTIONABLE_SIGNATURE_TSV = "actionableSignatures.tsv";
    private static final String ACTIONABLE_VARIANT_TSV = "actionableVariants.tsv";

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
        FusionEvidenceAnalyzer actionableFusion = FusionEvidenceAnalyzerFactory.loadFromActionableFusionTsv(basePath + ACTIONABLE_FUSION_TSV);
        GeneEvidenceAnalyzer actionableGene = GeneEvidenceAnalyzerFactory.loadFromActionableGeneTsv(basePath + ACTIONABLE_GENE_TSV);
        RangeEvidenceAnalyzer actionableRange = RangeEvidenceAnalyzerFactory.loadFromActionableRangeTsv(basePath + ACTIONABLE_RANGE_TSV);
        SignatureEvidenceAnalyzer actionableSignature =
                SignatureEvidenceAnalyzerFactory.loadFromActionableSignatureTsv(basePath + ACTIONABLE_SIGNATURE_TSV);
        VariantEvidenceAnalyzer actionableVariant = VariantEvidenceAnalyzerFactory.loadFromActionableVariantTsv(basePath + ACTIONABLE_VARIANT_TSV);

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

    @NotNull
    @VisibleForTesting
    public FusionEvidenceAnalyzer fusionEvidenceAnalyzer() {
        return fusionEvidenceAnalyzer;
    }

    @NotNull
    @VisibleForTesting
    public GeneEvidenceAnalyzer geneEvidenceAnalyzer() {
        return geneEvidenceAnalyzer;
    }

    @NotNull
    @VisibleForTesting
    RangeEvidenceAnalyzer rangeEvidenceAnalyzer() {
        return rangeEvidenceAnalyzer;
    }

    @NotNull
    @VisibleForTesting
    SignatureEvidenceAnalyzer signatureEvidenceAnalyzer() {
        return signatureEvidenceAnalyzer;
    }

    @NotNull
    @VisibleForTesting
    VariantEvidenceAnalyzer variantEvidenceAnalyzer() {
        return variantEvidenceAnalyzer;
    }
}
