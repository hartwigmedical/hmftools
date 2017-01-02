package com.hartwig.hmftools.patientreporter.slicing;

import org.jetbrains.annotations.NotNull;

public class GenomeRegion {

    @NotNull
    private final String chromosome;
    private final long start;
    private final long end;

    public GenomeRegion(@NotNull final String chromosome, final long start, final long end) {
        assert end >= start;

        this.chromosome = chromosome;
        this.start = start;
        this.end = end;
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

    long bases() {
        return 1 + end - start;
    }

    @Override
    public String toString() {
        return "GenomeRegion{" + "chromosome='" + chromosome + '\'' + ", start=" + start + ", end=" + end + '}';
    }
}
