package com.hartwig.hmftools.serve.sources.vicc.check;

import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import org.junit.Ignore;
import org.junit.Test;

public class GeneCheckerTest {

    @Test
    @Ignore
    public void canCorrectlyCheckGenes() {
        GeneChecker geneChecker = new GeneChecker();
        assertTrue(geneChecker.isValidGene("BRAF"));
        assertTrue(geneChecker.isValidGene("IGH"));

        assertFalse(geneChecker.isValidGene("I am not a gene"));
        assertFalse(geneChecker.isValidGene(null));
    }
}