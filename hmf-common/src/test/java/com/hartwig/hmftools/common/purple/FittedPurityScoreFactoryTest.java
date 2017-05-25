package com.hartwig.hmftools.common.purple;

import static com.hartwig.hmftools.common.purple.FittedPurityScoreFactory.isPolyclonal;

import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import org.junit.Test;

public class FittedPurityScoreFactoryTest
{

    @Test
    public void testPolyclonalProportion() {
        assertFalse(isPolyclonal(2.75));

        assertTrue(isPolyclonal(2.751));
        assertTrue(isPolyclonal(3.0));
        assertTrue(isPolyclonal(3.249));

        assertFalse(isPolyclonal(3.5));
    }

}
