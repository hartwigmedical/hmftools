package com.hartwig.hmftools.orange.algo.protect;

import org.jetbrains.annotations.NotNull;

public final class ProtectInterpretationTestFactory {

    private ProtectInterpretationTestFactory() {
    }

    @NotNull
    public static ImmutableProtectInterpretedData.Builder builder() {
        return ImmutableProtectInterpretedData.builder();
    }

    @NotNull
    public static ProtectInterpretedData createMinimalTestProtectData() {
        return builder().build();
    }
}
