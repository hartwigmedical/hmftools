package com.hartwig.hmftools.common.runcontext;

import org.jetbrains.annotations.NotNull;

public interface RunContext {

    @NotNull
    String runDirectory();

    @NotNull
    String setName();

    @NotNull
    String refSample();

    @NotNull
    String tumorSample();

    @NotNull
    String tumorBarcodeSample();
}
