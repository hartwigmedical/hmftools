package com.hartwig.hmftools.common.lims;

import static org.junit.Assert.assertEquals;

import org.junit.Test;

public class LimsCohortTest {

    @Test
    public void canExtractCohortTypeFromSample() {
        assertEquals(LimsCohort.CPCT, LimsCohort.fromCohort("CPCT", "CPCT01010001T"));
        assertEquals(LimsCohort.CPCT_PANCREAS, LimsCohort.fromCohort("CPCTPancreas", "CPCT01010000T"));
        assertEquals(LimsCohort.DRUP, LimsCohort.fromCohort("DRUP", "DRUP01010001T"));
        assertEquals(LimsCohort.DRUP_STAGE3, LimsCohort.fromCohort("DRUPstage3", "DRUP01010000T"));
        assertEquals(LimsCohort.CORE, LimsCohort.fromCohort("CORE", "CORE01010000T"));
        assertEquals(LimsCohort.CORELR02, LimsCohort.fromCohort("CORELR02", "CORELR020000T"));
        assertEquals(LimsCohort.CORERI02, LimsCohort.fromCohort("CORERI02", "CORERI020000T"));
        assertEquals(LimsCohort.CORELR11, LimsCohort.fromCohort("CORELR11", "CORELR110000T"));
        assertEquals(LimsCohort.CORESC11, LimsCohort.fromCohort("CORESC11", "CORESC110000T"));
        assertEquals(LimsCohort.WIDE, LimsCohort.fromCohort("WIDE", "WIDE010000T"));
        assertEquals(LimsCohort.COREDB, LimsCohort.fromCohort("COREDB", "COREDB010000T"));
        assertEquals(LimsCohort.COREDB, LimsCohort.fromCohort("COREDB", "WIDE01010000T"));
    }

    @Test(expected = IllegalStateException.class)
    public void hasUnknownCohortType() {
        LimsCohort.fromCohort("ABCD", "");
        LimsCohort.fromCohort("", "");
    }
}