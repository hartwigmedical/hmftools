package com.hartwig.hmftools.purple.config;

public class CommonConfig {
    private final String normalSample;
    private final String tumorSample;
    private final String outputDirectory;
    private final String runDirectory;
    private final String freecDirectory;
    private final boolean forceRecalculate;

    CommonConfig(final String normalSample, final String tumorSample, final String outputDirectory, final String runDirectory,
            final String freecDirectory, boolean forceRecalculate) {
        this.normalSample = normalSample;
        this.tumorSample = tumorSample;
        this.outputDirectory = outputDirectory;
        this.runDirectory = runDirectory;
        this.freecDirectory = freecDirectory;
        this.forceRecalculate = forceRecalculate;
    }

    public String refSample() {
        return normalSample;
    }

    public String tumorSample() {
        return tumorSample;
    }

    public String outputDirectory() {
        return outputDirectory;
    }

    public String runDirectory() {
        return runDirectory;
    }

    public String freecDirectory() {
        return freecDirectory;
    }

    public boolean forceRecalculate() {
        return forceRecalculate;
    }
}
