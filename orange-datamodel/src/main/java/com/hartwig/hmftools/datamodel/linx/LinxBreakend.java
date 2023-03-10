package com.hartwig.hmftools.datamodel.linx;

import com.hartwig.hmftools.datamodel.gene.TranscriptCodingType;
import com.hartwig.hmftools.datamodel.gene.TranscriptRegionType;
import com.hartwig.hmftools.datamodel.sv.LinxBreakendType;
import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;

@Value.Immutable
public abstract class LinxBreakend {
    public abstract int id();

    public abstract int svId();

    @NotNull
    public abstract String gene();

    @NotNull
    public abstract String transcriptId();

    public abstract boolean canonical();

    @NotNull
    public abstract String geneOrientation();

    public abstract boolean disruptive();

    public abstract boolean reportedDisruption();

    public abstract double undisruptedCopyNumber();

    @NotNull
    public abstract TranscriptRegionType regionType();

    @NotNull
    public abstract TranscriptCodingType codingType();

    public abstract int nextSpliceExonRank();

    // additional fields for patient report
    @NotNull
    public abstract LinxBreakendType type();

    @NotNull
    public abstract String chromosome();

    public abstract int orientation();

    public abstract int strand();

    @NotNull
    public abstract String chrBand();

    public abstract int exonUp();

    public abstract int exonDown();
    
    public abstract double junctionCopyNumber();
}
