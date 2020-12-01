package com.hartwig.hmftools.serve.sources.vicc.check;

import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import org.junit.Ignore;
import org.junit.Test;

public class GeneCheckerTest {

    @Test
    //TODO improve test
    public void canCorrectlyCheckGenes() {
        GeneChecker geneChecker = new GeneChecker();
        assertTrue(geneChecker.isValidGene("IGL", null, "Amplification"));
        assertTrue(geneChecker.isValidGene("IGH", null, "Deletion"));

        assertFalse(geneChecker.isValidGene("I am not a gene", null, "event"));
        assertFalse(geneChecker.isValidGene("gene", null, "event"));
    }
}