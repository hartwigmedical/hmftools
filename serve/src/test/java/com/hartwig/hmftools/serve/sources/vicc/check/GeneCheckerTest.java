package com.hartwig.hmftools.serve.sources.vicc.check;

import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import com.google.common.collect.Sets;

import org.junit.Test;

public class GeneCheckerTest {

    @Test
    public void canCorrectlyAssessGenes() {
        GeneChecker geneChecker = new GeneChecker(Sets.newHashSet("BRAF", "IGH"));

        assertTrue(geneChecker.isValidGene("BRAF", "Amplification"));
        assertTrue(geneChecker.isValidGene("IGH", "Amplification"));

        assertFalse(geneChecker.isValidGene("I am not a gene", "Amplification"));
        assertFalse(geneChecker.isValidGene(null, null));
    }
}