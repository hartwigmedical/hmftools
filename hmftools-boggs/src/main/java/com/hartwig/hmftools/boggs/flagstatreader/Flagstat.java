package com.hartwig.hmftools.boggs.flagstatreader;

import org.jetbrains.annotations.NotNull;

public class Flagstat {

    @NotNull
    private final String path;
    @NotNull
    private final Stats qcPassedReads;
    @NotNull
    private final Stats qcFailedReads;

    public Flagstat(@NotNull String path, @NotNull Stats qcPassedReads, @NotNull Stats qcFailedReads) {
        this.path = path;
        this.qcPassedReads = qcPassedReads;
        this.qcFailedReads = qcFailedReads;
    }

    @NotNull
    public Stats qcPassedReads() {
        return qcPassedReads;
    }

    @NotNull
    public Stats qcFailedReads() {
        return qcFailedReads;
    }

    @Override
    public String toString() {
        return "Flagstat{" +
                "path='" + path + '\'' +
                ", qcPassedReads=" + qcPassedReads +
                ", qcFailedReads=" + qcFailedReads +
                '}';
    }
}
