package com.hartwig.hmftools.boggs.flagstatreader;

import org.jetbrains.annotations.NotNull;

public class FlagStatData {

    @NotNull
    private final String path;
    @NotNull
    private final FlagStats qcPassedReads;
    @NotNull
    private final FlagStats qcFailedReads;

    public FlagStatData(@NotNull String path, @NotNull FlagStats qcPassedReads, @NotNull FlagStats qcFailedReads) {
        this.path = path;
        this.qcPassedReads = qcPassedReads;
        this.qcFailedReads = qcFailedReads;
    }

    @NotNull
    public String path() {
        return path;
    }

    @NotNull
    public FlagStats qcPassedReads() {
        return qcPassedReads;
    }

    @NotNull
    public FlagStats qcFailedReads() {
        return qcFailedReads;
    }

    @Override
    public String toString() {
        return "FlagStatData{" +
                "path='" + path + '\'' +
                ", qcPassedReads=" + qcPassedReads +
                ", qcFailedReads=" + qcFailedReads +
                '}';
    }
}
