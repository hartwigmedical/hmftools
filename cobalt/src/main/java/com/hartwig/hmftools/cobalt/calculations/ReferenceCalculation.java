package com.hartwig.hmftools.cobalt.calculations;

import static com.hartwig.hmftools.cobalt.CobaltConstants.ROLLING_MEDIAN_MAX_DISTANCE;
import static com.hartwig.hmftools.cobalt.CobaltConstants.ROLLING_MEDIAN_MIN_COVERAGE;

import java.util.List;

import com.hartwig.hmftools.cobalt.consolidation.ResultsConsolidator;
import com.hartwig.hmftools.cobalt.normalisers.DiploidNormaliser;
import com.hartwig.hmftools.cobalt.normalisers.ReadDepthStatisticsNormaliser;
import com.hartwig.hmftools.cobalt.normalisers.ResultsNormaliser;
import com.hartwig.hmftools.cobalt.targeted.CobaltScope;
import com.hartwig.hmftools.common.cobalt.MedianRatio;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion;

class ReferenceCalculation extends BamCalculation
{
    private final ResultsConsolidator mResultsConsolidator;

    ReferenceCalculation(final WindowStatuses mGenomeFilter, CobaltScope scope, RefGenomeVersion version,
            final ResultsConsolidator mResultsConsolidator)
    {
        super(mGenomeFilter, scope, version);
        this.mResultsConsolidator = mResultsConsolidator;
    }

    ReadDepthStatisticsNormaliser createReadDepthsNormaliser()
    {
        return Scope.medianByMeanNormaliser();
    }

    ResultsNormaliser createMegaBaseScaleNormaliser(RefGenomeVersion version)
    {
        return new DiploidNormaliser(ROLLING_MEDIAN_MAX_DISTANCE, ROLLING_MEDIAN_MIN_COVERAGE, version);
    }

    List<MedianRatio> medianRatios()
    {
        return ((DiploidNormaliser) MegaBaseScaleNormaliser).medianRatios();
    }

    @Override
    ResultsConsolidator consolidator()
    {
        // todo if null was passed in at construction, get the scope to build the consolidator
        return mResultsConsolidator;
    }
}
