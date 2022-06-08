package com.hartwig.hmftools.common.purple.interpretation;

import static org.junit.Assert.assertEquals;

import org.junit.Test;

public class CopyNumberInterpretationTest {

    @Test
    public void canRenderDisplay() {
        assertEquals("full gain", CopyNumberInterpretation.FULL_GAIN.display());
    }
}