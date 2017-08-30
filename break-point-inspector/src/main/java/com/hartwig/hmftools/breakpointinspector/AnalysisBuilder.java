package com.hartwig.hmftools.breakpointinspector;

import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SamReader;

class AnalysisBuilder {
    private SamReader refReader;
    private SAMFileWriter refWriter;
    private SamReader tumorReader;
    private SAMFileWriter tumorWriter;
    private int range = 500;
    private int[] extraUncertainty = { 0, 1, 5, 10, 20 };

    AnalysisBuilder setRefReader(final SamReader refReader) {
        this.refReader = refReader;
        return this;
    }

    AnalysisBuilder setRefWriter(final SAMFileWriter refWriter) {
        this.refWriter = refWriter;
        return this;
    }

    AnalysisBuilder setTumorReader(final SamReader tumorReader) {
        this.tumorReader = tumorReader;
        return this;
    }

    AnalysisBuilder setTumorWriter(final SAMFileWriter tumorWriter) {
        this.tumorWriter = tumorWriter;
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
        return new Analysis(refReader, refWriter, tumorReader, tumorWriter, range, extraUncertainty);
    }
}