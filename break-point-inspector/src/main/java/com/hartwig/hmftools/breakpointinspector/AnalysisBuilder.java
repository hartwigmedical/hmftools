package com.hartwig.hmftools.breakpointinspector;

import org.jetbrains.annotations.NotNull;

import htsjdk.samtools.SamReader;

class AnalysisBuilder {
    private SamReader refReader;
    private SamReader tumorReader;
    private int range = 500;
    private float contamination = 0;

    @NotNull
    AnalysisBuilder refReader(@NotNull SamReader refReader) {
        this.refReader = refReader;
        return this;
    }

    @NotNull
    AnalysisBuilder tumorReader(final SamReader tumorReader) {
        this.tumorReader = tumorReader;
        return this;
    }

    void range(final int range) {
        this.range = range;
    }

    void contaminationFraction(final float fraction) {
        this.contamination = fraction;
    }

    @NotNull
    Analysis createAnalysis() {
        return new Analysis(refReader, tumorReader, range, contamination);
    }
}