package com.hartwig.hmftools.cobalt.calculations;

import com.hartwig.hmftools.cobalt.targeted.CobaltScope;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion;

public class TumorBamCalculation extends BamCalculation
{
    public TumorBamCalculation(final GenomeFilter mGenomeFilter, CobaltScope scope, RefGenomeVersion version)
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
}
