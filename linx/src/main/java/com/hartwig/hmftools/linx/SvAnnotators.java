package com.hartwig.hmftools.linx;

import com.hartwig.hmftools.common.ensemblcache.EnsemblDataCache;
import com.hartwig.hmftools.linx.annotators.FragileSiteAnnotator;
import com.hartwig.hmftools.linx.annotators.IndelAnnotator;
import com.hartwig.hmftools.linx.annotators.KataegisAnnotator;
import com.hartwig.hmftools.linx.annotators.LineElementAnnotator;
import com.hartwig.hmftools.linx.annotators.PseudoGeneFinder;
import com.hartwig.hmftools.patientdb.dao.DatabaseAccess;

public class SvAnnotators
{
    public final FragileSiteAnnotator FragileSiteAnnotator;
    public final LineElementAnnotator LineElementAnnotator;
    public final KataegisAnnotator KataegisAnnotator;
    public final PseudoGeneFinder PseudoGeneFinder;
    public final IndelAnnotator IndelAnnotator;

    public SvAnnotators(
            final LinxConfig config, final EnsemblDataCache geneDataCache, final DatabaseAccess dbAccess,
            final CohortDataWriter cohortDataWriter)
    {
        FragileSiteAnnotator = new FragileSiteAnnotator();
        FragileSiteAnnotator.loadFragileSitesFile(config.FragileSiteFile);

        LineElementAnnotator = new LineElementAnnotator(config.ProximityDistance);
        LineElementAnnotator.loadLineElementsFile(config.LineElementFile);

        PseudoGeneFinder = new PseudoGeneFinder(geneDataCache);

        LineElementAnnotator.setPseudoGeneFinder(PseudoGeneFinder);

        KataegisAnnotator = new KataegisAnnotator(config.OutputDataPath);
        KataegisAnnotator.loadKataegisData(config.KataegisFile);

        IndelAnnotator = config.IndelAnnotation ? new IndelAnnotator(dbAccess, config) : null;
    }

    public void close()
    {
        KataegisAnnotator.close();

        if(IndelAnnotator != null)
            IndelAnnotator.close();
    }
}
