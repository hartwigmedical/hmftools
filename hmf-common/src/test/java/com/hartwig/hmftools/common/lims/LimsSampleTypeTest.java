package com.hartwig.hmftools.common.lims;

import static org.junit.Assert.assertEquals;

import org.junit.Test;

public class LimsSampleTypeTest {

    @Test
    public void canResolveLimsSampleTypes() {
        assertEquals(LimsStudy.CPCT, LimsStudy.fromSampleId("CPCT02990001T"));
        assertEquals(LimsStudy.WIDE, LimsStudy.fromSampleId("WIDEwideWIDE"));
        assertEquals(LimsStudy.CORE, LimsStudy.fromSampleId("CORE0299-CPCT01T"));
    }

    @Test
    public void unrecognizedSampleIdGivesOther() {
        assertEquals(LimsStudy.OTHER, LimsStudy.fromSampleId("Unknown"));
    }

    @Test
    public void lowerCaseGivesOther() {
        assertEquals(LimsStudy.OTHER, LimsStudy.fromSampleId("xxxdrupxxx"));
    }

    @Test
    public void doesNotStartWithGivesOther() {
        assertEquals(LimsStudy.OTHER, LimsStudy.fromSampleId("829COLO"));
    }
}