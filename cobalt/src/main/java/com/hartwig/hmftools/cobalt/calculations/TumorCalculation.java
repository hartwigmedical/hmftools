package com.hartwig.hmftools.cobalt.calculations;

import com.hartwig.hmftools.cobalt.consolidation.ResultsConsolidator;
import com.hartwig.hmftools.cobalt.normalisers.DoNothingNormaliser;
import com.hartwig.hmftools.cobalt.normalisers.ReadDepthStatisticsNormaliser;
import com.hartwig.hmftools.cobalt.normalisers.ResultsNormaliser;
import com.hartwig.hmftools.cobalt.targeted.CobaltScope;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion;

class TumorCalculation extends BamCalculation
{
    private ResultsConsolidator mResultsConsolidator;

    TumorCalculation(final WindowStatuses mGenomeFilter, CobaltScope scope)
    {
        super(mGenomeFilter, scope);
    }

    ReadDepthStatisticsNormaliser createReadDepthsNormaliser()
    {
        return Scope.medianByMeanNormaliser();
    }

    @Override
    ResultsNormaliser megaBaseScaleNormaliser()
    {
        return new DoNothingNormaliser();
    }

    @Override
    ResultsConsolidator consolidator()
    {
        if(mResultsConsolidator == null)
        {
            mResultsConsolidator = Scope.resultsConsolidator(MeanNormaliser.readDepthMedian());
        }
        return mResultsConsolidator;
    }
}
