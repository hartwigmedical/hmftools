package com.hartwig.hmftools.patientreporter.slicing;

public class GenomeRegion {

    private final long start;
    private final long end;

    public GenomeRegion(final long start, final long end) {
        assert end >= start;

        this.start = start;
        this.end = end;
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
        return "GenomeRegion{" + "start=" + start + ", end=" + end + '}';
    }
}
