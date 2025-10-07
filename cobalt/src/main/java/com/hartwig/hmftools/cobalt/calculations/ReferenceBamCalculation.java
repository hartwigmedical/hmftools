package com.hartwig.hmftools.cobalt.calculations;

import static com.hartwig.hmftools.cobalt.CobaltConstants.ROLLING_MEDIAN_MAX_DISTANCE;
import static com.hartwig.hmftools.cobalt.CobaltConstants.ROLLING_MEDIAN_MIN_COVERAGE;

import java.util.List;

import com.hartwig.hmftools.cobalt.targeted.CobaltScope;
import com.hartwig.hmftools.common.cobalt.MedianRatio;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion;

public class ReferenceBamCalculation extends BamCalculation
{
    public ReferenceBamCalculation(final GenomeFilter mGenomeFilter, CobaltScope scope, RefGenomeVersion version)
    {
        super(mGenomeFilter, scope, version);
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
}
