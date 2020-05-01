package com.hartwig.hmftools.common.lims;

import static org.junit.Assert.assertEquals;

import org.junit.Test;

public class LimsCoreCohortTest {

    @Test
    public void canResolveLimsCoreCohorts() {
        assertEquals(LimsCoreCohort.CORELR02, LimsCoreCohort.fromSampleId("CORELR020000T"));
        assertEquals(LimsCoreCohort.CORELR11, LimsCoreCohort.fromSampleId("CORELR110000T"));
    }

    @Test
    public void unrecognizedSampleIdGivesNonCore() {
        assertEquals(LimsCoreCohort.NON_CORE, LimsCoreCohort.fromSampleId("Unknown"));
    }

    @Test
    public void lowerCaseGivesNonCore() {
        assertEquals(LimsCoreCohort.NON_CORE, LimsCoreCohort.fromSampleId("xxxcorexxx"));
    }
}