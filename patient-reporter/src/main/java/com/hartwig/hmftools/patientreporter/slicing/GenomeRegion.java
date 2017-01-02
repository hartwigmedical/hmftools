package com.hartwig.hmftools.patientreporter.slicing;

import org.jetbrains.annotations.NotNull;

class GenomeRegion implements Comparable<GenomeRegion> {

    private final long start;
    private final long end;

    GenomeRegion(final long start, final long end) {
        assert end >= start;

        this.start = start;
        this.end = end;
    }

    long start() {
        return start;
    }

    long end() {
        return end;
    }

    long bases() {
        return 1 + end - start;
    }

    @Override
    public String toString() {
        return "GenomeRegion{" + "start=" + start + ", end=" + end + '}';
    }

    @Override
    public int compareTo(@NotNull final GenomeRegion other) {
        if (start() < other.start()) {
            return -1;
        } else if (start() == other.start()) {
            return 0;
        }
        return 1;
    }
}
