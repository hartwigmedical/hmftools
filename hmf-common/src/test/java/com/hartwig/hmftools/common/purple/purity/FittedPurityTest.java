package com.hartwig.hmftools.common.purple.purity;

import static org.junit.Assert.assertEquals;

import org.jetbrains.annotations.NotNull;
import org.junit.Test;

public class FittedPurityTest {

    private static final double EPSILON = 1e-10;

    @Test
    public void testPloidy() {
        assertEquals(7, create(0.10, 0.80).ploidy(), EPSILON);
    }

    @NotNull
    private static FittedPurity create(final double purity, final double normFactor) {
        return ImmutableFittedPurity.builder().score(0).diploidProportion(0).modelBAFDeviation(0).purity(
                purity).normFactor(normFactor).build();
    }
}
