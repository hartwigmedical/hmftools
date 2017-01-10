package com.hartwig.hmftools.healthchecker.context;

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

    boolean isSomaticRun();
}
