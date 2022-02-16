package com.hartwig.hmftools.serve.refgenome;

import java.io.IOException;

import com.google.common.io.Resources;
import com.hartwig.hmftools.common.fusion.KnownFusionCache;

import org.jetbrains.annotations.NotNull;


public class RefGenomeManagerFactoryTest  {

    private static final String FUSION_DATA_37 = Resources.getResource("known_fusion_data/known_fusion_data.37.csv").getPath();

    private RefGenomeManagerFactoryTest() {
    }

    @NotNull
    public static KnownFusionCache knownFusionCache() {
        try {
            return RefGenomeManagerFactory.buildKnownFusionCacheFromFile(FUSION_DATA_37);
        } catch (IOException e) {
            throw new IllegalStateException("Could not load test fusion data cache");
        }
    }
}