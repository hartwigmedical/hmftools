package com.hartwig.hmftools.common.variant.structural.linx;

import org.immutables.value.Value;

@Value.Immutable
public abstract class LinxBreakend
{
    public abstract int id();
    public abstract int svId();
    public abstract boolean isStart();
    public abstract String gene();
    public abstract String transcriptId();
    public abstract boolean canonical();
    public abstract boolean isUpstream();
    public abstract boolean disruptive();
    public abstract boolean reportedDisruption();
    public abstract double undisruptedCopyNumber();
    public abstract String regionType();
    public abstract String codingContext();
    public abstract String biotype();
    public abstract int exonBasePhase();
    public abstract int nextSpliceExonRank();
    public abstract int nextSpliceExonPhase();
    public abstract int nextSpliceDistance();
    public abstract int totalExonCount();

}
