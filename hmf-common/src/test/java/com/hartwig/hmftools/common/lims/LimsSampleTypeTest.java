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

    @Test (expected = IllegalStateException.class)
    public void unrecognizedSampleIdGivesException() {
        LimsSampleType.fromSampleId("Unknown");
    }

    @Test (expected = IllegalStateException.class)
    public void lowerCaseGivesException() {
        LimsSampleType.fromSampleId("xxxdrupxxx");
    }

    @Test (expected = IllegalStateException.class)
    public void doesNotStartWithGivesException() {
        LimsSampleType.fromSampleId("829COLO");
    }
}