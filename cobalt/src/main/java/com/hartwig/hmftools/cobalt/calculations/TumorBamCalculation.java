package com.hartwig.hmftools.cobalt.calculations;

import com.hartwig.hmftools.cobalt.targeted.CobaltScope;

public class TumorBamCalculation extends BamCalculation
{
    public TumorBamCalculation(final GenomeFilter mGenomeFilter, CobaltScope scope)
    {
        super(mGenomeFilter, scope);
    }

    ResultsNormaliser finalMeanNormaliser()
    {
        return mScope.finalNormaliser();
    }

    ResultsNormaliser diploidNormaliser()
    {
        return new DoNothingNormaliser();
    }
}
