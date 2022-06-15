package com.hartwig.hmftools.orange.algo.linx;

import org.jetbrains.annotations.NotNull;

public final class LinxInterpreterTestFactory {

    private LinxInterpreterTestFactory() {
    }

    @NotNull
    public static LinxInterpretedData createMinimalTestLinxData() {
        return ImmutableLinxInterpretedData.builder().build();
    }
}
