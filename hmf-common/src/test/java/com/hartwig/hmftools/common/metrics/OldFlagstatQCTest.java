package com.hartwig.hmftools.common.metrics;

import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import org.jetbrains.annotations.NotNull;
import org.junit.Test;

public class OldFlagstatQCTest
{
    @Test
    public void canEvaluateFlagstatQC()
    {
        OldFlagstat qcPass = withMappedProportion(0.99);
        OldFlagstat qcFail = withMappedProportion(0.01);

        assertTrue(FlagstatQC.pass(qcPass));
        assertFalse(FlagstatQC.pass(qcFail));
    }

    @NotNull
    private static OldFlagstat withMappedProportion(double mappedProportion)
    {
        return ImmutableOldFlagstat.builder().from(FlagstatTestFactory.createMinimalTestFlagstat()).mappedProportion(mappedProportion).build();
    }
}