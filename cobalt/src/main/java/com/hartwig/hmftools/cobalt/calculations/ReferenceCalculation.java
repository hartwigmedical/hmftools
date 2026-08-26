package com.hartwig.hmftools.cobalt.calculations;

import static com.hartwig.hmftools.cobalt.CobaltConstants.ROLLING_MEDIAN_MAX_DISTANCE;
import static com.hartwig.hmftools.cobalt.CobaltConstants.ROLLING_MEDIAN_MIN_COVERAGE;

import java.util.Collections;
import java.util.List;

import com.hartwig.hmftools.cobalt.consolidation.ResultsConsolidator;
import com.hartwig.hmftools.cobalt.normalisers.DiploidNormaliser;
import com.hartwig.hmftools.cobalt.normalisers.DoNothingNormaliser;
import com.hartwig.hmftools.cobalt.normalisers.ReadDepthStatisticsNormaliser;
import com.hartwig.hmftools.cobalt.normalisers.ResultsNormaliser;
import com.hartwig.hmftools.cobalt.targeted.CobaltScope;
import com.hartwig.hmftools.common.cobalt.MedianRatio;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion;

class ReferenceCalculation extends BamCalculation
{
    private final ResultsConsolidator mResultsConsolidator;
    private final ResultsNormaliser mMegaBaseDiploidNormaliser;

    ReferenceCalculation(
            final WindowStatuses mGenomeFilter, final CobaltScope scope, final RefGenomeVersion version,
            final ResultsConsolidator resultsConsolidator, boolean requireMegaBaseScaleDiploidNorm)
    {
        super(mGenomeFilter, scope);
        mResultsConsolidator = resultsConsolidator;

        if(requireMegaBaseScaleDiploidNorm)
        {
            mMegaBaseDiploidNormaliser = new DiploidNormaliser(ROLLING_MEDIAN_MAX_DISTANCE, ROLLING_MEDIAN_MIN_COVERAGE, version);
        }
        else
        {
            mMegaBaseDiploidNormaliser = new DoNothingNormaliser();
        }
    }

    ReadDepthStatisticsNormaliser createReadDepthsNormaliser()
    {
        return Scope.medianByMeanNormaliser();
    }

    List<MedianRatio> medianRatios()
    {
        if(mMegaBaseDiploidNormaliser instanceof DiploidNormaliser)
            return ((DiploidNormaliser)mMegaBaseDiploidNormaliser).medianRatios();
        else
            return Collections.emptyList();
    }

    @Override
    public ResultsNormaliser megaBaseScaleNormaliser() { return mMegaBaseDiploidNormaliser; }

    @Override
    ResultsConsolidator consolidator()
    {
        // todo if null was passed in at construction, get the scope to build the consolidator
        return mResultsConsolidator;
    }
}
