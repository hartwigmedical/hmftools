package com.hartwig.hmftools.datamodel.linx;

import com.hartwig.hmftools.datamodel.gene.TranscriptCodingType;
import com.hartwig.hmftools.datamodel.gene.TranscriptRegionType;
import com.hartwig.hmftools.datamodel.sv.LinxBreakendType;
import org.immutables.value.Value;

@Value.Immutable
public abstract class LinxBreakend {
    public abstract int svId();
    public abstract String gene();
    public abstract String transcriptId();
    public abstract boolean canonical();
    public abstract String geneOrientation();
    public abstract boolean disruptive();
    public abstract boolean reportedDisruption();
    public abstract double undisruptedCopyNumber();
    public abstract TranscriptRegionType regionType();
    public abstract TranscriptCodingType codingType();
    public abstract int nextSpliceExonRank();

    // additional fields for patient report
    public abstract LinxBreakendType type();
    public abstract String chromosome();
    public abstract int orientation();
    public abstract int strand();
    public abstract String chrBand();
    public abstract int exonUp();
    public abstract int exonDown();
    public abstract double junctionCopyNumber();
}
