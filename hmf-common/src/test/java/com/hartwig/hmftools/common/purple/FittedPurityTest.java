package com.hartwig.hmftools.common.purple;

import static org.junit.Assert.assertEquals;

import org.junit.Test;

public class FittedPurityTest {

    private static final double EPSILON = 1e-10;

    @Test
    public void testPloidy() {
        assertEquals(7, create(0.10, 0.80).ploidy(), EPSILON);
    }


    private FittedPurity create(double purity, double normFactor) {
        return ImmutableFittedPurity.builder()
                .score(0)
                .diplodProportion(0)
                .modelBAFDeviation(0)
                .purity(purity)
                .normFactor(normFactor)
                .build();
    }

}
