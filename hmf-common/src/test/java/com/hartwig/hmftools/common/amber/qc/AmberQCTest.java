package com.hartwig.hmftools.common.amber.qc;

import static com.hartwig.hmftools.common.amber.qc.AmberQCStatus.FAIL;
import static com.hartwig.hmftools.common.amber.qc.AmberQCStatus.PASS;
import static com.hartwig.hmftools.common.amber.qc.AmberQCStatus.WARN;

import static org.junit.Assert.assertEquals;

import org.jetbrains.annotations.NotNull;
import org.junit.Test;

public class AmberQCTest {

    @Test
    public void testStatus() {
        assertEquals(FAIL, create(0.511).status());
        assertEquals(WARN, create(0.510).status());
        assertEquals(WARN, create(0.501).status());
        assertEquals(PASS, create(0.500).status());
        assertEquals(PASS, create(0.487).status());
        assertEquals(WARN, create(0.486).status());
        assertEquals(WARN, create(0.480).status());
        assertEquals(FAIL, create(0.479).status());
    }


    @NotNull
    private static AmberQC create(double meanBAF) {
        return ImmutableAmberQC.builder().meanBAF(meanBAF).build();
    }

}
