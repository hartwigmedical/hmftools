package com.hartwig.hmftools.protect.actionability_v2;

import java.io.File;
import java.io.IOException;
import java.util.List;

import com.hartwig.hmftools.protect.actionability_v2.fusion.ActionableFusion;
import com.hartwig.hmftools.protect.actionability_v2.gene.ActionableGene;

import org.jetbrains.annotations.NotNull;

public class ActionabilityAnalyzer {

    private static final String ACTIONABLE_FUSION_TSV = "actionableFusion.tsv";
    private static final String ACTIONABLE_GENE_TSV = "actionableGene.tsv";
    private static final String ACTIONABLE_RANGE_TSV = "actionableRange.tsv";
    private static final String ACTIONABLE_SIGNATURE_TSV = "actionableSignature.tsv";
    private static final String ACTIONABLE_VARIANT_TSV = "actionableVariant.tsv";

    @NotNull
    public static void fromKnowledgebase(@NotNull String knowledgebaseDirectory) throws
            IOException {

        String basePath = knowledgebaseDirectory + File.separator;
        List<ActionableFusion> actionableFusion = ReadActionabilityFiles.loadFromFileFusion(basePath + ACTIONABLE_FUSION_TSV);
        List<ActionableGene> actionableGene = ReadActionabilityFiles.loadFromFileGene(basePath + ACTIONABLE_GENE_TSV);

    }
}
