package com.hartwig.hmftools.common.lims;

import static org.junit.Assert.*;

import org.junit.Test;

public class LimsCohortTest {

    @Test
    public void canExtractCohortTypeFromSample() {
        assertEquals(LimsCohort.CPCT, LimsCohort.fromCohort("CPCT", "CPCT010000001T"));
        assertEquals(LimsCohort.CPCT, LimsCohort.fromCohort("CPCTpancreas", "CPCT010000001T"));
        assertEquals(LimsCohort.DRUP, LimsCohort.fromCohort("DRUP", "DRUP10000001T"));
        assertEquals(LimsCohort.DRUP, LimsCohort.fromCohort("DRUPstage3", "DRUP010000001T"));
        assertEquals(LimsCohort.CORE, LimsCohort.fromCohort("CORE", "CORE010000001T"));
        assertEquals(LimsCohort.CORELR02, LimsCohort.fromCohort("CORELR02", "CORELR020000T"));
        assertEquals(LimsCohort.CORERI02, LimsCohort.fromCohort("CORERI02", "CORERI020000T"));
        assertEquals(LimsCohort.CORELR11, LimsCohort.fromCohort("CORELR11", "CORELR110000T"));
        assertEquals(LimsCohort.CORESC11, LimsCohort.fromCohort("CORESC11", "CORESC110000T"));
        assertEquals(LimsCohort.WIDE, LimsCohort.fromCohort("WIDE", "WIDE10000001T"));
        assertEquals(LimsCohort.COREDB01, LimsCohort.fromCohort("COREDB01", "COREDB010000T"));
    }

    @Test(expected = IllegalStateException.class)
    public void hasUnknownCohortType() {

        LimsCohort.fromCohort("ABCD", "ACDV01010000T");
        LimsCohort.fromCohort("", "ACDV01010000T");
    }
}