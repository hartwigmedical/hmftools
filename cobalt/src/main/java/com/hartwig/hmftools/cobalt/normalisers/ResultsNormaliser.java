package com.hartwig.hmftools.cobalt.normalisers;

import com.hartwig.hmftools.cobalt.calculations.BamRatio;

public interface ResultsNormaliser
{
    void recordValue(BamRatio bamRatio);

    default void dataCollectionFinished()
    {
    }

    void normalise(BamRatio bamRatio);
}
