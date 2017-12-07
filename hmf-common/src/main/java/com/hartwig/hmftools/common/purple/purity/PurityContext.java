package com.hartwig.hmftools.common.purple.purity;

import com.hartwig.hmftools.common.purple.gender.Gender;

import org.immutables.value.Value;

@Value.Immutable
public abstract class PurityContext {

    public abstract String version();
    
    public abstract Gender gender();

    public abstract FittedPurity bestFit();

    public abstract FittedPurityStatus status();

    public abstract FittedPurityScore score();

    public abstract double polyClonalProportion();
}
