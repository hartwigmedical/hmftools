package com.hartwig.hmftools.amber.contamination;

public record VafConsistencyCheckResult(double gini, int totalWeightInBand, int totalWeightAcrossAllVafValues)
{
}
