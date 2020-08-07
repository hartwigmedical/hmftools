package com.hartwig.hmftools.common.sigs;

import org.immutables.value.Value;

@Value.Immutable
public abstract class SignatureAllocation
{
    public static final String SIG_UNALLOCATED = "UNALLOC";
    public static final String SIG_EXCESS = "EXCESS";

    public abstract String signature();
    public abstract double allocation();
    public abstract double percent();
}
