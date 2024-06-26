package com.hartwig.hmftools.linx;

import com.hartwig.hmftools.common.ensemblcache.EnsemblDataCache;
import com.hartwig.hmftools.linx.annotators.FragileSiteAnnotator;
import com.hartwig.hmftools.linx.annotators.LineElementAnnotator;
import com.hartwig.hmftools.linx.annotators.PseudoGeneFinder;

public class SvAnnotators
{
    public final FragileSiteAnnotator FragileSiteAnnotator;
    public final LineElementAnnotator LineElementAnnotator;
    public final PseudoGeneFinder PseudoGeneFinder;

    public SvAnnotators(final LinxConfig config, final EnsemblDataCache geneDataCache)
    {
        FragileSiteAnnotator = new FragileSiteAnnotator();
        FragileSiteAnnotator.loadFragileSitesFile(config.FragileSiteFile);

        LineElementAnnotator = new LineElementAnnotator(config.ProximityDistance);
        LineElementAnnotator.loadLineElementsFile(config.LineElementFile);

        PseudoGeneFinder = new PseudoGeneFinder(geneDataCache);

        LineElementAnnotator.setPseudoGeneFinder(PseudoGeneFinder);
    }
}
