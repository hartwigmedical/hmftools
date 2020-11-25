package com.hartwig.hmftools.common.serve.classification;

import static org.junit.Assert.assertEquals;

import java.util.Map;

import com.google.common.collect.Maps;

import org.jetbrains.annotations.NotNull;
import org.junit.Test;

public class EventClassifierTest {

    private static final String HOTSPOT = "hotspot";
    private static final String FUSION_PAIR = "fusion pair";
    private static final String MULTIPLE = "multiple";

    @Test
    public void canDetermineMutationTypes() {
        EventClassifier classifier = new EventClassifier(buildTestMatcherMap());

        assertEquals(MutationType.HOTSPOT, classifier.determineType("any", HOTSPOT));
        assertEquals(MutationType.FUSION_PAIR, classifier.determineType("any", FUSION_PAIR));

        // Events with multiple types should be UNKNOWN.
        assertEquals(MutationType.UNKNOWN, classifier.determineType("any", MULTIPLE));

        assertEquals(MutationType.UNKNOWN, classifier.determineType("any", "any"));
    }

    @NotNull
    private static Map<MutationType, EventMatcher> buildTestMatcherMap() {
        Map<MutationType, EventMatcher> map = Maps.newHashMap();
        map.put(MutationType.HOTSPOT, (gene, event) -> event.equals(HOTSPOT));
        map.put(MutationType.FUSION_PAIR, (gene, event) -> event.equals(FUSION_PAIR));
        map.put(MutationType.COMPLEX, (gene, event) -> event.equals(MULTIPLE));
        map.put(MutationType.COMBINED, (gene, event) -> event.equals(MULTIPLE));
        return map;
    }
}