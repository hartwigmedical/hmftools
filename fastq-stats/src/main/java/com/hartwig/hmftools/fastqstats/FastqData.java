package com.hartwig.hmftools.fastqstats;

import org.jetbrains.annotations.NotNull;

class FastqData {

    private final long yield;
    private final long q30;

    FastqData(final long yield, final long q30) {
        this.q30 = q30;
        this.yield = yield;
    }

    long yield() {
        return yield;
    }

    long q30() {
        return q30;
    }

    double q30Percentage() {
        return q30 * 100.0 / yield;
    }

    @NotNull
    FastqData add(@NotNull final FastqData other) {
        return new FastqData(yield + other.yield, q30 + other.q30);
    }
}
