package com.hartwig.hmftools.common.variant.structural.linx;

import org.immutables.value.Value;

@Value.Immutable
public abstract class LinxFusion
{
    public abstract int fivePrimeBreakendId();
    public abstract int threePrimeBreakendId();
    public abstract String name();
    public abstract boolean reported();
    public abstract String reportedType();
    public abstract boolean phased();
    public abstract int chainLength();
    public abstract int chainLinks();
    public abstract boolean chainTerminated();
    public abstract String domainsKept();
    public abstract String domainsLost();
    public abstract int skippedExonsUp();
    public abstract int skippedExonsDown();
    public abstract int fusedExonUp();
    public abstract int fusedExonDown();

}
