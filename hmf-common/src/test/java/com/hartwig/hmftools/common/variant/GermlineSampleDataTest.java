package com.hartwig.hmftools.common.variant;

import static org.junit.Assert.assertEquals;

import org.apache.logging.log4j.util.Strings;
import org.junit.Test;

public class GermlineSampleDataTest {

    private static final double EPSILON = 1.0e-10;

    @Test
    public void canCalculateAlleleFrequency() {
        final GermlineSampleData data = new GermlineSampleData(Strings.EMPTY, 20, 10);
        assertEquals(0.5, data.alleleFrequency(), EPSILON);
    }
}