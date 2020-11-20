package com.hartwig.hmftools.vicc.annotation;

import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.serve.classification.EventMatcher;

import org.junit.Test;

public class DeletionClassifierTest {

    @Test
    public void canAssessWhetherEventIsDeletion() {
        EventMatcher classifier = DeletionClassifier.create(Lists.newArrayList());

        assertTrue(classifier.matches("CDKN2A", "CDKN2A del"));
        assertTrue(classifier.matches("CDKN2A", "CDKN2A dec exp"));

        assertFalse(classifier.matches("BRAF", "V600E"));
    }
}