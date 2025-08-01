package com.hartwig.hmftools.common.utils;

import java.io.File;
import java.net.URISyntaxException;
import java.net.URL;
import java.util.Objects;

import com.hartwig.hmftools.common.ensemblcache.EnsemblDataCache;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion;

public class EnsemblMini
{
    public static File ensemblMiniDataDir()
    {
        URL resourceUrl = EnsemblMini.class.getClassLoader().getResource("ensembl_mini");
        Objects.requireNonNull(resourceUrl, "Resource 'ensembl_mini' not found in classpath");
        try
        {
            return new File(resourceUrl.toURI());
        }
        catch(URISyntaxException e)
        {
            throw new RuntimeException(e);
        }
    }

    public static EnsemblDataCache ensemblMiniDataCache()
    {
        EnsemblDataCache cache = new EnsemblDataCache(ensemblMiniDataDir().getAbsolutePath(), RefGenomeVersion.V38);
        cache.setRequiredData(true, true, true, false);
        cache.load(false);
        cache.createTranscriptIdMap();
        return cache;
    }
}
