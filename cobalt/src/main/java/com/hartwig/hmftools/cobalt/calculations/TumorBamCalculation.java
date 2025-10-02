package com.hartwig.hmftools.cobalt.calculations;

import com.hartwig.hmftools.cobalt.targeted.CobaltScope;

public class TumorBamCalculation extends BamCalculation
{
    public TumorBamCalculation(final GenomeFilter mGenomeFilter, CobaltScope scope)
    {
        super(mGenomeFilter, scope);
    }

    ReadDepthStatisticsNormaliser createReadDepthsNormaliser()
    {
        return mScope.medianByMeanNormaliser();
    }

    ResultsNormaliser createMegaBaseScaleNormaliser()
    {
        return new DoNothingNormaliser();
    }

    @Override
    ResultsNormaliser createFinalNormaliser()
    {
        return mScope.finalNormaliser();
    }
}
