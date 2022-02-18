package com.hartwig.hmftools.serve;

import java.io.IOException;

import com.google.common.io.Resources;
import com.hartwig.hmftools.common.ensemblcache.EnsemblDataCache;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion;
import com.hartwig.hmftools.serve.refgenome.EnsemblDataCacheLoader;

import org.jetbrains.annotations.NotNull;

public final class EnsemblDataCacheTestFactory {

    private static final String ENSEMBL_DATA_DIR_37 = Resources.getResource("ensembl_data_cache/v37").getPath();

    private EnsemblDataCacheTestFactory() {
    }

    @NotNull
    public static EnsemblDataCache create37() {
        try {
            return EnsemblDataCacheLoader.load(ENSEMBL_DATA_DIR_37, RefGenomeVersion.V37);
        } catch (IOException e) {
            throw new IllegalStateException("Could not load test ensembl cache");
        }
    }
}