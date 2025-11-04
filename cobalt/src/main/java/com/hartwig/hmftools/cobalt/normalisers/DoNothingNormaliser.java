package com.hartwig.hmftools.cobalt.normalisers;

import com.hartwig.hmftools.cobalt.calculations.BamRatio;

public class DoNothingNormaliser implements ResultsNormaliser
{
    @Override
    public void recordValue(final BamRatio bamRatio)
    {
    }

    @Override
    public void normalise(BamRatio bamRatio)
    {
    }
}
