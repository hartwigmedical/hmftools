package com.hartwig.hmftools.amber.contamination;

import java.util.List;

public record VafConsistencyCheckResult<T extends Comparable<T>>(double unevenDistributionCost, int totalWeightInBand,
                                                                 int totalWeightAcrossAllVafValues,
                                                                 List<CategoryEvidence<T>> categoryEvidence)
{
    double proportion()
    {
        return totalWeightInBand / (double) totalWeightAcrossAllVafValues;
    }
}
