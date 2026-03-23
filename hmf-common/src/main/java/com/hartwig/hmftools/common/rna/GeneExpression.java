package com.hartwig.hmftools.common.rna;

import com.hartwig.hmftools.common.purple.ReportedStatus;

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
    public abstract ReportedStatus reportedStatus();

    public static final double NO_CANCER_AVAILABLE_VALUE = -1;

}
