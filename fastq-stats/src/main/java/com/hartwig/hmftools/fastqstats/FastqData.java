package com.hartwig.hmftools.fastqstats;

import org.jetbrains.annotations.NotNull;

public class FastqData {
    private final long yield;
    private final long q30;

    public FastqData(long yield, long q30) {
        this.q30 = q30;
        this.yield = yield;
    }

    public long getYield() {
        return yield;
    }

    public long getQ30() {
        return q30;
    }

    @NotNull
    public FastqData add(@NotNull FastqData other) {
        return new FastqData(yield + other.yield, q30 + other.q30);
    }
}
