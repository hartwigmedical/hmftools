package com.hartwig.hmftools.patientreporter.slicing;

import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public class GenomeRegion implements Comparable<GenomeRegion> {

    @NotNull
    private final String chromosome;
    private final long start;
    private final long end;
    @Nullable
    private final String annotation;

    public GenomeRegion(@NotNull final String chromosome, final long start, final long end) {
        this(chromosome, start, end, null);
    }

    GenomeRegion(@NotNull final String chromosome, final long start, final long end,
            @Nullable final String annotation) {
        assert end >= start;

        this.chromosome = chromosome;
        this.start = start;
        this.end = end;
        this.annotation = annotation;
    }

    @NotNull
    public String chromosome() {
        return chromosome;
    }

    @Nullable
    public String annotation() {
        return annotation;
    }

    public long start() {
        return start;
    }

    public long end() {
        return end;
    }

    public long bases() {
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
