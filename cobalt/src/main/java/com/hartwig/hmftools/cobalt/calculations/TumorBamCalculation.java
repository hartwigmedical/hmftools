package com.hartwig.hmftools.cobalt.calculations;

import com.hartwig.hmftools.cobalt.targeted.TargetRegions;

public class TumorBamCalculation extends BamCalculation
{
    public TumorBamCalculation(final GenomeFilter mGenomeFilter, TargetRegions targetRegions)
    {
        super(mGenomeFilter, targetRegions);
    }

    ResultsNormaliser finalMeanNormaliser()
    {
        return new UnityNormaliser();
    }

    ResultsNormaliser diploidNormaliser()
    {
        return new DoNothingNormaliser();
    }
}
