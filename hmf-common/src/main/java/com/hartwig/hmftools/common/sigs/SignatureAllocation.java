package com.hartwig.hmftools.common.sigs;

import org.immutables.value.Value;

@Value.Immutable
public abstract class SignatureAllocation
{
    public abstract String signature();
    public abstract double allocation();
    public abstract double percent();
}
