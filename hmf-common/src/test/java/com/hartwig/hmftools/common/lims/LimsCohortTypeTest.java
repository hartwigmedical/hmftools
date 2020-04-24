package com.hartwig.hmftools.common.lims;

import static org.junit.Assert.*;

import org.junit.Test;

public class LimsCohortTypeTest {

    @Test
    public void canResolveLimsCohortTypes() {
        assertEquals(LimsCohortType.CORELR02, LimsCohortType.fromSampleId("CORELR020000T"));
        assertEquals(LimsCohortType.CORELR11, LimsCohortType.fromSampleId("CORELR110000T"));
    }

    @Test
    public void unrecognizedSampleIdGivesOther() {
        assertEquals(LimsCohortType.OTHER, LimsCohortType.fromSampleId("Unknown"));
    }

    @Test
    public void lowerCaseGivesOther() {
        assertEquals(LimsCohortType.OTHER, LimsCohortType.fromSampleId("xxxcorexxx"));
    }

}