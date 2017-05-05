package com.hartwig.hmftools.common.variant;

import org.junit.Test;

import static org.junit.Assert.assertEquals;

public class GermlineSampleDataTest {

    private static final double EPSILON = 1.0e-10;

    @Test
    public void canCalculateAlleleFrequency() {
        final GermlineSampleData data = ImmutableGermlineSampleData.of("1/1", 20, 10, 35);

        assertEquals(0.5, data.alleleFrequency(), EPSILON);
    }
}