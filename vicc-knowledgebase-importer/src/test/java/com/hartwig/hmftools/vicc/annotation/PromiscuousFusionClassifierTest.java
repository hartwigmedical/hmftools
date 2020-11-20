package com.hartwig.hmftools.vicc.annotation;

import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.serve.classification.EventClassifier;

import org.junit.Test;

public class PromiscuousFusionClassifierTest {

    @Test
    public void canAssessWhetherEventIsPromiscuousFusion() {
        EventClassifier classifier = PromiscuousFusionClassifier.create(Lists.newLinkedList());

        assertTrue(classifier.matches("BRAF", "BRAF fusion"));
        assertTrue(classifier.matches("ROS1", "ROS1 REARRANGEMENT"));

        assertFalse(classifier.matches("BRAF", "V600E"));
    }
}