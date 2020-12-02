package com.hartwig.hmftools.serve.fusion;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.serve.Knowledgebase;

import org.jetbrains.annotations.NotNull;
import org.junit.Test;

public class FusionFunctionsTest {

    @Test
    public void canConsolidateFusionPairs() {
        String gene1 = "gene1";
        String gene2 = "gene2";
        KnownFusionPair fusion1 = ImmutableKnownFusionPair.builder()
                .geneUp(gene1)
                .geneDown(gene2)
                .addSources(Knowledgebase.VICC_ONCOKB, Knowledgebase.VICC_CIVIC)
                .build();
        KnownFusionPair fusion2 = ImmutableKnownFusionPair.builder().geneUp(gene1).geneDown(gene2).addSources(Knowledgebase.VICC_CGI).build();
        KnownFusionPair fusion3 = ImmutableKnownFusionPair.builder().geneUp(gene2).geneDown(gene1).addSources(Knowledgebase.VICC_CGI).build();

        List<KnownFusionPair> consolidated = FusionFunctions.consolidate(Lists.newArrayList(fusion1, fusion2, fusion3));
        assertEquals(2, consolidated.size());

        KnownFusionPair gene1Fusion = findByGeneUp(consolidated, gene1);
        assertEquals(3, gene1Fusion.sources().size());
        assertTrue(gene1Fusion.sources().contains(Knowledgebase.VICC_CGI));
        assertTrue(gene1Fusion.sources().contains(Knowledgebase.VICC_ONCOKB));
        assertTrue(gene1Fusion.sources().contains(Knowledgebase.VICC_CIVIC));

        KnownFusionPair gene2Fusion = findByGeneUp(consolidated, gene2);
        assertEquals(1, gene2Fusion.sources().size());
        assertTrue(gene2Fusion.sources().contains(Knowledgebase.VICC_CGI));
    }

    @NotNull
    private static KnownFusionPair findByGeneUp(@NotNull List<KnownFusionPair> fusionPairs, @NotNull String gene) {
        for (KnownFusionPair fusionPair : fusionPairs) {
            if (fusionPair.geneUp().equals(gene)) {
                return fusionPair;
            }
        }

        throw new IllegalStateException("Could not find gene in fusion pairs: " + gene);
    }
}