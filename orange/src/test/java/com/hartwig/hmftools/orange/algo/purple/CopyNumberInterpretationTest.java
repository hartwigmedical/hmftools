package com.hartwig.hmftools.orange.algo.purple;

import static org.junit.Assert.assertEquals;

import org.junit.Test;

public class CopyNumberInterpretationTest {

    @Test
    public void canRenderDisplay() {
        assertEquals("full gain", CopyNumberInterpretation.FULL_GAIN.display());
    }
}