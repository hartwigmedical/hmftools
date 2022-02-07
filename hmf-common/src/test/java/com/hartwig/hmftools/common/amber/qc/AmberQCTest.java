package com.hartwig.hmftools.common.amber.qc;

import static com.hartwig.hmftools.common.amber.qc.AmberQCStatus.FAIL;
import static com.hartwig.hmftools.common.amber.qc.AmberQCStatus.PASS;
import static com.hartwig.hmftools.common.amber.qc.AmberQCStatus.WARN;

import static org.junit.Assert.assertEquals;

import org.jetbrains.annotations.NotNull;
import org.junit.Test;

public class AmberQCTest {

    @Test
    public void testMeanStatus() {
        assertEquals(FAIL, createWithMean(0.511).status());
        assertEquals(WARN, createWithMean(0.510).status());
        assertEquals(WARN, createWithMean(0.501).status());
        assertEquals(PASS, createWithMean(0.500).status());
        assertEquals(PASS, createWithMean(0.487).status());
        assertEquals(WARN, createWithMean(0.486).status());
        assertEquals(WARN, createWithMean(0.480).status());
        assertEquals(FAIL, createWithMean(0.479).status());
    }

    @Test
    public void testContaminationStatus() {
        assertEquals(PASS, createWithContamination(0).status());
        assertEquals(WARN, createWithContamination(0.001).status());
        assertEquals(WARN, createWithContamination(0.10).status());
        assertEquals(FAIL, createWithContamination(0.101).status());
        assertEquals(FAIL, createWithContamination(1).status());
    }

    @NotNull
    private static AmberQC createWithMean(double meanBAF) {
        return ImmutableAmberQC.builder()
                .meanBAF(meanBAF)
                .contamination(0)
                .consanguinityProportion(0)
                .uniparentalDisomy(null).build();
    }

    @NotNull
    private static AmberQC createWithContamination(double contamination) {
        return ImmutableAmberQC.builder()
                .meanBAF(0.5)
                .contamination(contamination)
                .consanguinityProportion(0)
                .uniparentalDisomy(null).build();
    }
}
