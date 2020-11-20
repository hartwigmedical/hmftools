package com.hartwig.hmftools.vicc.annotation;

import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.serve.classification.EventMatcher;

import org.junit.Test;

public class GeneRangeExonClassifierTest {

    @Test
    public void canAssessWhetherEventIsGeneRangeExon() {
        EventMatcher classifier = GeneRangeExonClassifier.create(Lists.newArrayList());

        assertTrue(classifier.matches("EGFR", "RARE EX 18-21 MUT"));
        assertTrue(classifier.matches("EGFR", "EGFR exon 19 deletions"));
        assertTrue(classifier.matches("AXSL1", "EXON 12 MUTATION"));

        assertFalse(classifier.matches("BRAF", "V600E"));
    }
}