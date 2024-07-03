package com.hartwig.hmftools.purple.fitting;

import static com.hartwig.hmftools.purple.fitting.FittedPurityScoreFactory.isPolyclonal;

import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import org.junit.Test;

public class FittedPurityScoreFactoryTest
{

    @Test
    public void testPolyclonalProportion()
    {
        assertTrue(isPolyclonal(2.749));

        assertFalse(isPolyclonal(2.75));
        assertFalse(isPolyclonal(2.751));
        assertFalse(isPolyclonal(3.0));
        assertFalse(isPolyclonal(3.249));
        assertFalse(isPolyclonal(3.25));

        assertTrue(isPolyclonal(3.251));
    }
}
