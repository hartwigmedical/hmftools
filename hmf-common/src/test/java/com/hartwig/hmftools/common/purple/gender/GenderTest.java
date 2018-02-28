package com.hartwig.hmftools.common.purple.gender;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNull;

import com.google.common.collect.ArrayListMultimap;
import com.google.common.collect.Multimap;

import org.junit.Test;

public class GenderTest {

    @Test
    public void testGender() {

        assertEquals(Gender.MALE,Gender.fromRatio(createRatios(), x -> x));
        assertEquals(Gender.MALE,Gender.fromRatio(createRatios(-1), x -> x));
        assertEquals(Gender.MALE, Gender.fromRatio(createRatios(0.75), x -> x));
        assertEquals(Gender.FEMALE, Gender.fromRatio(createRatios(0.76), x -> x));
    }

    private static Multimap<String, Double> createRatios(double... ratios) {
        final Multimap<String, Double> map = ArrayListMultimap.create();
        for (double ratio : ratios) {
            map.put("X", ratio);
        }

        return map;
    }

}
