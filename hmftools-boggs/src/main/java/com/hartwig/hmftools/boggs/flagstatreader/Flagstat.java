package com.hartwig.hmftools.boggs.flagstatreader;

import org.jetbrains.annotations.NotNull;

public class Flagstat {

    @NotNull
    private final Stats qcPassedReads;
    @NotNull
    private final Stats qcFailedReads;

    Flagstat(@NotNull Stats qcPassedReads, @NotNull Stats qcFailedReads) {
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
}
