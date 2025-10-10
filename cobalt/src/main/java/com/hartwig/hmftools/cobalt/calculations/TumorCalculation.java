package com.hartwig.hmftools.cobalt.calculations;

import com.hartwig.hmftools.cobalt.consolidation.ResultsConsolidator;
import com.hartwig.hmftools.cobalt.targeted.CobaltScope;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion;

public class TumorCalculation extends BamCalculation
{
    private ResultsConsolidator mResultsConsolidator;
    public TumorCalculation(final GenomeFilter mGenomeFilter, CobaltScope scope, RefGenomeVersion version)
    {
        super(mGenomeFilter, scope, version);
    }

    ReadDepthStatisticsNormaliser createReadDepthsNormaliser()
    {
        return Scope.medianByMeanNormaliser();
    }

    ResultsNormaliser createMegaBaseScaleNormaliser(RefGenomeVersion version)
    {
        return new DoNothingNormaliser();
    }

    @Override
    ResultsConsolidator consolidator()
    {
        if (mResultsConsolidator == null)
        {
            mResultsConsolidator = Scope.resultsConsolidator(MeanNormaliser.readDepthMedian());
        }
        return mResultsConsolidator;
    }
}
