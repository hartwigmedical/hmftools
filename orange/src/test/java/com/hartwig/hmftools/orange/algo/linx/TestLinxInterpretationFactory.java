package com.hartwig.hmftools.orange.algo.linx;

import org.jetbrains.annotations.NotNull;

public final class TestLinxInterpretationFactory {

    private TestLinxInterpretationFactory() {
    }

    @NotNull
    public static LinxInterpretedData createMinimalTestLinxData() {
        return ImmutableLinxInterpretedData.builder().build();
    }
}
