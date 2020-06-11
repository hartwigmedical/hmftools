package com.hartwig.hmftools.common.variant.structural.linx;

import org.immutables.value.Value;

@Value.Immutable
public abstract class LinxLink
{
    public abstract int clusterId();
    public abstract int chainId();
    public abstract String chainIndex();
    public abstract int chainCount();
    public abstract int lowerSvId();
    public abstract int upperSvId();
    public abstract boolean lowerBreakendIsStart();
    public abstract boolean upperBreakendIsStart();
    public abstract String chromosome();
    public abstract String arm();
    public abstract boolean assembled();
    public abstract int traversedSVCount();
    public abstract long length();
    public abstract double junctionCopyNumber();
    public abstract double junctionCopyNumberUncertainty();
    public abstract String pseudogeneInfo();
}
