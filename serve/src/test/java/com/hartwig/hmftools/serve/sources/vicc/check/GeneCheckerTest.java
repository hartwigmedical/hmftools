package com.hartwig.hmftools.serve.sources.vicc.check;

import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import org.junit.Ignore;
import org.junit.Test;

public class GeneCheckerTest {

    @Test
    //TODO improve test
    @Ignore
    public void canCorrectlyCheckGenes() {
        GeneChecker geneChecker = new GeneChecker();
        assertTrue(geneChecker.isValidGene("BRAF", null, "Amplification", null));
        assertTrue(geneChecker.isValidGene("IGH", null, "Deletion", null));

        assertFalse(geneChecker.isValidGene("I am not a gene", null, "event", null));
        assertFalse(geneChecker.isValidGene("gene", null, "event", null));
    }
}