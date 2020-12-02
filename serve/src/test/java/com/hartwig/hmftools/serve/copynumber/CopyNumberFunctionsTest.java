package com.hartwig.hmftools.serve.copynumber;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.serve.Knowledgebase;

import org.jetbrains.annotations.NotNull;
import org.junit.Test;

public class CopyNumberFunctionsTest {

    @Test
    public void canConsolidateCopyNumbers() {
        String gene1 = "gene1";
        String gene2 = "gene2";
        KnownCopyNumber copyNumber1 = ampBuilder().gene(gene1).addSources(Knowledgebase.VICC_ONCOKB, Knowledgebase.VICC_CIVIC).build();
        KnownCopyNumber copyNumber2 = ampBuilder().gene(gene1).addSources(Knowledgebase.VICC_CGI).build();
        KnownCopyNumber copyNumber3 = ampBuilder().gene(gene2).addSources(Knowledgebase.VICC_CGI).build();

        List<KnownCopyNumber> consolidated = CopyNumberFunctions.consolidate(Lists.newArrayList(copyNumber1, copyNumber2, copyNumber3));
        assertEquals(2, consolidated.size());

        KnownCopyNumber gene1CopyNumber = findByGene(consolidated, gene1);
        assertEquals(3, gene1CopyNumber.sources().size());
        assertTrue(gene1CopyNumber.sources().contains(Knowledgebase.VICC_CGI));
        assertTrue(gene1CopyNumber.sources().contains(Knowledgebase.VICC_ONCOKB));
        assertTrue(gene1CopyNumber.sources().contains(Knowledgebase.VICC_CIVIC));

        KnownCopyNumber gene2CopyNumber = findByGene(consolidated, gene2);
        assertEquals(1, gene2CopyNumber.sources().size());
        assertTrue(gene2CopyNumber.sources().contains(Knowledgebase.VICC_CGI));
    }

    @NotNull
    private static ImmutableKnownCopyNumber.Builder ampBuilder() {
        return ImmutableKnownCopyNumber.builder().type(CopyNumberType.AMPLIFICATION);
    }

    @NotNull
    private static KnownCopyNumber findByGene(@NotNull List<KnownCopyNumber> copyNumbers, @NotNull String gene) {
        for (KnownCopyNumber copyNumber : copyNumbers) {
            if (copyNumber.gene().equals(gene)) {
                return copyNumber;
            }
        }

        throw new IllegalStateException("Could not find gene in copy numbers: " + gene);
    }
}