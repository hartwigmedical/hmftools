package com.hartwig.hmftools.common.lims;

import static org.junit.Assert.*;

import org.junit.Test;

public class LimsCohortTypeTest {

    @Test
    public void canResolveLimsCohortTypes() {
        assertEquals(LimsCoreCohort.CORELR02, LimsCoreCohort.fromSampleId("CORELR020000T"));
        assertEquals(LimsCoreCohort.CORELR11, LimsCoreCohort.fromSampleId("CORELR110000T"));
    }

    @Test
    public void unrecognizedSampleIdGivesOther() {
        assertEquals(LimsCoreCohort.OTHER, LimsCoreCohort.fromSampleId("Unknown"));
    }

    @Test
    public void lowerCaseGivesOther() {
        assertEquals(LimsCoreCohort.OTHER, LimsCoreCohort.fromSampleId("xxxcorexxx"));
    }

}