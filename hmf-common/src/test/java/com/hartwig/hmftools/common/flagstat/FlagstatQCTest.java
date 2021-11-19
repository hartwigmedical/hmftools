package com.hartwig.hmftools.common.flagstat;

import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import org.jetbrains.annotations.NotNull;
import org.junit.Test;

public class FlagstatQCTest {

    @Test
    public void canEvaluateFlagstatQC() {
        Flagstat qcPass = withMappedProportion(0.99);
        Flagstat qcFail = withMappedProportion(0.01);

        assertTrue(FlagstatQC.pass(qcPass));
        assertFalse(FlagstatQC.pass(qcFail));
    }

    @NotNull
    private static Flagstat withMappedProportion(double mappedProportion) {
        return ImmutableFlagstat.builder().from(FlagstatTestFactory.createMinimalTestFlagstat()).mappedProportion(mappedProportion).build();
    }
}