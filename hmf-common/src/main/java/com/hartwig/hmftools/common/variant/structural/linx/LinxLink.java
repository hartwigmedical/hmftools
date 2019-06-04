package com.hartwig.hmftools.common.variant.structural.linx;

import org.immutables.value.Value;

@Value.Immutable
public abstract class LinxLink
{
    public abstract int clusterId();
    public abstract int chainId();
    public abstract int chainIndex();
    public abstract int chainCount();
    public abstract int lowerBreakendId();
    public abstract int upperBreakendId();
    public abstract boolean lowerBreakendIsStart();
    public abstract boolean upperBreakendIsStart();
    public abstract String arm();
    public abstract boolean assembled();
    public abstract int traversedSVCount();
    public abstract long length();
    public abstract double ploidy();
    public abstract String pseudogeneInfo();
}
