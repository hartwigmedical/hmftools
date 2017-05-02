package com.hartwig.hmftools.common.variant;

import static org.junit.Assert.assertEquals;

import org.apache.logging.log4j.util.Strings;
import org.junit.Test;

public class GermlineSampleDataTest {

    private static final double EPSILON = 1.0e-10;
    public static ImmutableGermlineSampleData.Builder GERMLINE_SAMPLE_BUILDER = ImmutableGermlineSampleData.builder()
            .genoType("1/1")
            .alleleReadCount(10)
            .totalReadCount(20)
            .combinedDepth(0);

    @Test
    public void canCalculateAlleleFrequency() {
        final GermlineSampleData data = GERMLINE_SAMPLE_BUILDER
                .alleleReadCount(10)
                .totalReadCount(20)
                .build();

        assertEquals(0.5, data.alleleFrequency(), EPSILON);
    }
}