package com.hartwig.hmftools.cobalt.calculations;

public interface ResultsNormaliser
{
    void recordValue(BamRatio bamRatio);

    default void recordsAllAdded()
    {
    }

    void applyNormalisation(BamRatio bamRatio);
}
