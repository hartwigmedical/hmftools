package com.hartwig.hmftools.common.serve.classification.matchers;

import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import java.util.Set;

import com.google.common.collect.Sets;

import org.junit.Test;

public class CharacteristicMatcherTest {

    private static final Set<String> CHARACTERISTIC_EVENTS = Sets.newHashSet("Characteristic");

    @Test
    public void canAssessWhetherEventIsCharacteristic() {
        EventMatcher matcher = new CharacteristicMatcher(CHARACTERISTIC_EVENTS);

        assertTrue(matcher.matches("-", "Characteristic"));
        assertTrue(matcher.matches("-", "This is a Characteristic event"));

        assertFalse(matcher.matches("BRAF", "V600E"));
    }
}