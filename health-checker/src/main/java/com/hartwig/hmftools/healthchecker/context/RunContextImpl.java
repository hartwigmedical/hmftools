package com.hartwig.hmftools.healthchecker.context;

import org.jetbrains.annotations.NotNull;

class RunContextImpl implements RunContext {

    @NotNull
    private final String runDirectory;
    @NotNull
    private final String runName;
    @NotNull
    private final String refSample;
    @NotNull
    private final String tumorSample;
    private final boolean hasPassedTests;
    private final boolean isSomaticRun;

    RunContextImpl(@NotNull final String runDirectory, @NotNull final String runName, @NotNull final String refSample,
            @NotNull final String tumorSample, final boolean hasPassedTests, final boolean isSomaticRun) {
        this.runDirectory = runDirectory;
        this.runName = runName;
        this.refSample = refSample;
        this.tumorSample = tumorSample;
        this.hasPassedTests = hasPassedTests;
        this.isSomaticRun = isSomaticRun;
    }

    @NotNull
    @Override
    public String runDirectory() {
        return runDirectory;
    }

    @NotNull
    @Override
    public String runName() {
        return runName;
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

    @Override
    public boolean hasPassedTests() {
        return hasPassedTests;
    }

    @Override
    public boolean isSomaticRun() {
        return isSomaticRun;
    }
}

