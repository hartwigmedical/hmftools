package com.hartwig.hmftools.common.amber.qc;

import static com.hartwig.hmftools.common.amber.qc.AmberQCStatus.FAIL;
import static com.hartwig.hmftools.common.amber.qc.AmberQCStatus.PASS;
import static com.hartwig.hmftools.common.amber.qc.AmberQCStatus.WARN;

import static org.junit.Assert.assertEquals;

import org.jetbrains.annotations.NotNull;
import org.junit.Test;

public class AmberQCTest {

    @Test
    public void testContaminationStatus() {
        assertEquals(PASS, createWithContamination(0).status());
        assertEquals(WARN, createWithContamination(0.001).status());
        assertEquals(WARN, createWithContamination(0.10).status());
        assertEquals(FAIL, createWithContamination(0.101).status());
        assertEquals(FAIL, createWithContamination(1).status());
    }

    @NotNull
    private static AmberQC createWithContamination(double contamination) {
        return ImmutableAmberQC.builder()
                .contamination(contamination)
                .consanguinityProportion(0)
                .uniparentalDisomy(null).build();
    }
}
