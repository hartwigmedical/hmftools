package com.hartwig.hmftools.common.copynumber;

import org.jetbrains.annotations.NotNull;

public class CopyNumber {

    private static final int EXPECTED_VALUE = 2;

    @NotNull
    private final String chromosome;
    private final long start;
    private final long end;
    private final int value;

    public CopyNumber(@NotNull final String chromosome, final long start, final long end, final int value) {
        assert end >= start;
        assert value >= 0;
        assert value != EXPECTED_VALUE;

        this.chromosome = chromosome;
        this.start = start;
        this.end = end;
        this.value = value;
    }

    @NotNull
    public String chromosome() {
        return chromosome;
    }

    public long start() {
        return start;
    }

    public long end() {
        return end;
    }

    public long basesAffected() {
        return 1 + end - start;
    }

    public int value() {
        return value;
    }

    public boolean isGain() {
        return value > EXPECTED_VALUE;
    }

    public boolean isLoss() {
        return value < EXPECTED_VALUE;
    }
}
