package com.hartwig.hmftools.common.context;

import org.jetbrains.annotations.NotNull;

class RunContextImpl implements RunContext {

    @NotNull
    private final String runDirectory;
    @NotNull
    private final String setName;
    @NotNull
    private final String refSample;
    @NotNull
    private final String tumorSample;

    RunContextImpl(@NotNull final String runDirectory, @NotNull final String setName, @NotNull final String refSample,
            @NotNull final String tumorSample) {
        this.runDirectory = runDirectory;
        this.setName = setName;
        this.refSample = refSample;
        this.tumorSample = tumorSample;
    }

    @NotNull
    @Override
    public String runDirectory() {
        return runDirectory;
    }

    @NotNull
    @Override
    public String setName() {
        return setName;
    }

    @NotNull
    @Override
    public String refSample() {
        return refSample;
    }

    @NotNull
    @Override
    public String tumorSample() {
        return tumorSample;
    }
}

