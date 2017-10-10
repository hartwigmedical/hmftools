package com.hartwig.hmftools.common.purple.purity;

import org.immutables.value.Value;

@Value.Immutable
public abstract class SomaticPeak {

    public abstract double adjustedVAF();

    public abstract int count();
}
