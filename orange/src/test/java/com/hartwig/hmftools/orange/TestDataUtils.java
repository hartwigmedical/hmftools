package com.hartwig.hmftools.orange;

import com.google.common.io.Resources;
import com.hartwig.hmftools.common.ensemblcache.EnsemblDataCache;
import com.hartwig.hmftools.common.genome.chromosome.CytoBands;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion;

public final class TestDataUtils
{
    private static final String ENSEMBL_DATA_CACHE_PATH = Resources.getResource("ensembl").getPath();
    public static final CytoBands CYTO_BANDS = new CytoBands(RefGenomeVersion.V37);

    public static EnsemblDataCache createDummyCache()
    {
        EnsemblDataCache cache = new EnsemblDataCache(ENSEMBL_DATA_CACHE_PATH, RefGenomeVersion.V37);
        return cache;
    }

    public static EnsemblDataCache loadTestCache()
    {
        EnsemblDataCache cache = new EnsemblDataCache(ENSEMBL_DATA_CACHE_PATH, RefGenomeVersion.V37);
        cache.setRequireNonEnsemblTranscripts();
        cache.load(false);
        return cache;
    }
}
