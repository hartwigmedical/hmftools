package com.hartwig.hmftools.serve.extraction.util;

import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.serve.classification.EventType;

import org.junit.Test;

public class GeneCheckerTest {

    @Test
    public void canCorrectlyAssessGenes() {
        GeneChecker geneChecker = new GeneChecker(Sets.newHashSet("BRAF", "IGH"));

        assertTrue(geneChecker.isValidGene("BRAF", EventType.AMPLIFICATION));
        assertTrue(geneChecker.isValidGene("IGH", EventType.PROMISCUOUS_FUSION));

        assertFalse(geneChecker.isValidGene("I am not a gene", EventType.UNKNOWN));
        assertFalse(geneChecker.isValidGene(null, EventType.UNKNOWN));
    }
}