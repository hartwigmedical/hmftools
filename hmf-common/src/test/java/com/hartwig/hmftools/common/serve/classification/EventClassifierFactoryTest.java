package com.hartwig.hmftools.common.serve.classification;

import static org.junit.Assert.assertEquals;

import org.junit.Test;

public class EventClassifierFactoryTest {

    @Test
    public void canBuildEventClassifier() {
        EventClassifier classifier = EventClassifierFactory.buildClassifier(ImmutableEventClassifierConfig.builder()
                .proteinAnnotationExtractor(event -> event)
                .build());

        assertEquals(MutationType.UNKNOWN, classifier.determineType("gene", "mutation"));
    }
}