package com.hartwig.hmftools.common.purple.purity;

import org.immutables.value.Value;

@Value.Immutable
public abstract class FittedPurityScore {

    public abstract double minPurity();

    public abstract double maxPurity();

    public abstract double minPloidy();

    public abstract double maxPloidy();

    public abstract double minDiploidProportion();

    public abstract double maxDiploidProportion();
}
