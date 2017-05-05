package com.hartwig.hmftools.common.convoy;

import org.immutables.value.Value;

@Value.Immutable
public abstract class ConvoyPurity {

    public abstract double purity();

    public abstract double normFactor();

    public abstract double score();

    public abstract double modelBAFDeviation();

    public abstract double diplodProportion();
}
