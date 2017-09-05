package com.hartwig.hmftools.breakpointinspector;

import htsjdk.samtools.SamReader;

class AnalysisBuilder {
    private SamReader refReader;
    private SamReader tumorReader;
    private int range = 500;
    private int[] extraUncertainty = { 1, 5, 10, 20 };

    AnalysisBuilder setRefReader(final SamReader refReader) {
        this.refReader = refReader;
        return this;
    }

    AnalysisBuilder setTumorReader(final SamReader tumorReader) {
        this.tumorReader = tumorReader;
        return this;
    }

    AnalysisBuilder setRange(final int range) {
        this.range = range;
        return this;
    }

    AnalysisBuilder setExtraUncertainty(final int[] extraUncertainty) {
        this.extraUncertainty = extraUncertainty;
        return this;
    }

    Analysis createAnalysis() {
        return new Analysis(refReader, tumorReader, range, extraUncertainty);
    }
}