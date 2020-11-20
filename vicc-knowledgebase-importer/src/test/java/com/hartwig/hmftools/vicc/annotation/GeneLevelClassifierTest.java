package com.hartwig.hmftools.vicc.annotation;

import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.serve.classification.EventClassifier;

import org.junit.Test;

public class GeneLevelClassifierTest {

    @Test
    public void canAssessWhetherEventIsGeneLevelEvent() {
        EventClassifier classifier = GeneLevelClassifier.create(Lists.newArrayList());

        assertTrue(classifier.matches("AKT1", "AKT1 act mut"));
        assertTrue(classifier.matches("AKT1", "AKT1"));
        assertTrue(classifier.matches("ATK1", "positive"));
        assertTrue(classifier.matches("TP53", "biallelic inactivation"));

        assertFalse(classifier.matches("NPM", "Exon 12 mutation"));
        assertFalse(classifier.matches("BRAF", "V600E"));
    }
}