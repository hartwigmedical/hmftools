package com.hartwig.hmftools.cobalt.calculations;

public interface ResultsNormaliser
{
    void recordValue(BamRatio bamRatio);

    default void dataCollectionFinished()
    {
    }

    void normalise(BamRatio bamRatio);
}
