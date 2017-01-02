package com.hartwig.hmftools.common.io.dir;

import org.jetbrains.annotations.NotNull;

public interface RunContext {

    @NotNull
    String runDirectory();

    @NotNull
    String runName();

    @NotNull
    String refSample();

    @NotNull
    String tumorSample();

    boolean hasPassedTests();
}
