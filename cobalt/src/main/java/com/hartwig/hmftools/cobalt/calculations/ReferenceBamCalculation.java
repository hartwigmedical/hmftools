package com.hartwig.hmftools.cobalt.calculations;

import static com.hartwig.hmftools.cobalt.CobaltConstants.ROLLING_MEDIAN_MAX_DISTANCE;
import static com.hartwig.hmftools.cobalt.CobaltConstants.ROLLING_MEDIAN_MIN_COVERAGE;

import java.util.List;

import com.hartwig.hmftools.cobalt.targeted.CobaltScope;
import com.hartwig.hmftools.common.cobalt.MedianRatio;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion;

public class ReferenceBamCalculation extends BamCalculation
{
    private final DiploidNormaliser mDiploidNormaliser;

    public ReferenceBamCalculation(final GenomeFilter mGenomeFilter, CobaltScope scope, RefGenomeVersion version)
    {
        super(mGenomeFilter, scope);
        mDiploidNormaliser = new DiploidNormaliser(ROLLING_MEDIAN_MAX_DISTANCE, ROLLING_MEDIAN_MIN_COVERAGE, version);
    }

    ReadDepthStatisticsNormaliser createReadDepthsNormaliser()
    {
        return mScope.medianByMeanNormaliser();
    }

    ResultsNormaliser createMegaBaseScaleNormaliser()
    {
        return mDiploidNormaliser;
    }

    @Override
    ResultsNormaliser createFinalNormaliser()
    {
        return mScope.finalNormaliser();
    }

    List<MedianRatio> medianRatios()
    {
        return mDiploidNormaliser.medianRatios();
    }
}
