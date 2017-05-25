package com.hartwig.hmftools.common.purple;

import org.immutables.value.Value;

@Value.Immutable
public abstract class FittedPurityScore {

    public abstract double polyclonalProportion();

    public abstract double minPurity();

    public abstract double maxPurity();

    public abstract double minPloidy();

    public abstract double maxPloidy();
}
