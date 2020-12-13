package com.hartwig.hmftools.common.lims;

import static org.junit.Assert.assertEquals;

import org.junit.Test;

public class LimsCohortTest {

    @Test
    public void canExtractCohortTypeFromSample() {
        assertEquals(LimsCohort.CPCT, LimsCohort.fromLimsCohortString("CPCT", "CPCT01010001T"));
        assertEquals(LimsCohort.CPCT_PANCREAS, LimsCohort.fromLimsCohortString("CPCTPancreas", "CPCT01010000T"));
        assertEquals(LimsCohort.DRUP, LimsCohort.fromLimsCohortString("DRUP", "DRUP01010001T"));
        assertEquals(LimsCohort.DRUP_STAGE3, LimsCohort.fromLimsCohortString("DRUPstage3", "DRUP01010000T"));
        assertEquals(LimsCohort.CORE, LimsCohort.fromLimsCohortString("CORE", "CORE01010000T"));
        assertEquals(LimsCohort.CORELR02, LimsCohort.fromLimsCohortString("CORELR02", "CORELR020000T"));
        assertEquals(LimsCohort.CORERI02, LimsCohort.fromLimsCohortString("CORERI02", "CORERI020000T"));
        assertEquals(LimsCohort.CORELR11, LimsCohort.fromLimsCohortString("CORELR11", "CORELR110000T"));
        assertEquals(LimsCohort.CORESC11, LimsCohort.fromLimsCohortString("CORESC11", "CORESC110000T"));
        assertEquals(LimsCohort.WIDE, LimsCohort.fromLimsCohortString("WIDE", "WIDE010000T"));
        assertEquals(LimsCohort.COREDB, LimsCohort.fromLimsCohortString("COREDB", "COREDB010000T"));
        assertEquals(LimsCohort.COREDB, LimsCohort.fromLimsCohortString("COREDB", "WIDE01010000T"));
    }

    @Test(expected = IllegalStateException.class)
    public void hasUnknownCohortType() {
        LimsCohort.fromLimsCohortString("ABCD", "");
        LimsCohort.fromLimsCohortString("", "");
    }
}