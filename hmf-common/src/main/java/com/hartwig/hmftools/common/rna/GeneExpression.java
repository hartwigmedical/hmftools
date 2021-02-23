package com.hartwig.hmftools.common.rna;

import org.immutables.value.Value;

@Value.Immutable
public abstract class GeneExpression
{
    public abstract String geneName();
    public abstract double tpm();
    public abstract int splicedFragments();
    public abstract int unsplicedFragments();
    public abstract double medianTpmCancer();
    public abstract double percentileCancer();
    public abstract double medianTpmCohort();
    public abstract double percentileCohort();
}
