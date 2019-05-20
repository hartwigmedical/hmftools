package com.hartwig.hmftools.common.lims;

import static org.junit.Assert.assertEquals;

import org.junit.Test;

public class LimsSampleTypeTest {

    @Test
    public void canResolveLimsSampleTypes() {
        assertEquals(LimsSampleType.CPCT, LimsSampleType.fromSampleId("CPCT02990001T"));
        assertEquals(LimsSampleType.WIDE, LimsSampleType.fromSampleId("WIDEwideWIDE"));
        assertEquals(LimsSampleType.CORE, LimsSampleType.fromSampleId("CORE0299-CPCT01T"));
    }

    @Test
    public void unrecognizedSampleIdGivesOther() {
        assertEquals(LimsSampleType.OTHER, LimsSampleType.fromSampleId("Unknown"));
    }

    @Test
    public void lowerCaseGivesOther() {
        assertEquals(LimsSampleType.OTHER, LimsSampleType.fromSampleId("xxxdrupxxx"));
    }

    @Test
    public void doesNotStartWithGivesOther() {
        assertEquals(LimsSampleType.OTHER, LimsSampleType.fromSampleId("829COLO"));
    }
}