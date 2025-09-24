package com.hartwig.hmftools.cobalt.calculations;

import static com.hartwig.hmftools.cobalt.CobaltConstants.ROLLING_MEDIAN_MAX_DISTANCE;
import static com.hartwig.hmftools.cobalt.CobaltConstants.ROLLING_MEDIAN_MIN_COVERAGE;

import com.hartwig.hmftools.cobalt.targeted.TargetRegions;

public class ReferenceBamCalculation extends BamCalculation
{
    private final DiploidNormaliser mDiploidNormaliser = new DiploidNormaliser(ROLLING_MEDIAN_MAX_DISTANCE, ROLLING_MEDIAN_MIN_COVERAGE);

    public ReferenceBamCalculation(final GenomeFilter mGenomeFilter, TargetRegions targetRegions)
    {
        super(mGenomeFilter, targetRegions);
    }

    ResultsNormaliser finalMeanNormaliser()
    {
        return new DoNothingNormaliser();
    }

    ResultsNormaliser diploidNormaliser()
    {
        return mDiploidNormaliser;
    }

}
