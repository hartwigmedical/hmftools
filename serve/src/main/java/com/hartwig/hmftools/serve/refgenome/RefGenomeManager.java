package com.hartwig.hmftools.serve.refgenome;

import java.util.Map;

import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.genome.genepanel.HmfGenePanelSupplier;
import com.hartwig.hmftools.common.serve.RefGenomeVersion;
import com.hartwig.hmftools.serve.ServeConfig;

import org.jetbrains.annotations.NotNull;

public class RefGenomeManager {

    @NotNull
    private final Map<RefGenomeVersion, RefGenomeResource> refGenomeResourceMap;

    @NotNull
    public static RefGenomeManager buildFromServeConfig(@NotNull ServeConfig serveConfig) {
        RefGenomeResource refGenome37 = ImmutableRefGenomeResource.builder()
                .fastaFile(serveConfig.refGenome37FastaFile())
                .canonicalTranscriptPerGeneMap(HmfGenePanelSupplier.allGenesMap37())
                .putChainToOtherRefGenomeMap(RefGenomeVersion.V38, serveConfig.refGenome37To38Chain())
                .build();

        RefGenomeResource refGenome38 = ImmutableRefGenomeResource.builder()
                .fastaFile(serveConfig.refGenome38FastaFile())
                .canonicalTranscriptPerGeneMap(HmfGenePanelSupplier.allGenesMap38())
                .putChainToOtherRefGenomeMap(RefGenomeVersion.V37, serveConfig.refGenome38To37Chain())
                .build();

        Map<RefGenomeVersion, RefGenomeResource> refGenomeResourceMap = Maps.newHashMap();
        refGenomeResourceMap.put(RefGenomeVersion.V37, refGenome37);
        refGenomeResourceMap.put(RefGenomeVersion.V38, refGenome38);
        return new RefGenomeManager(refGenomeResourceMap);
    }

    private RefGenomeManager(@NotNull final Map<RefGenomeVersion, RefGenomeResource> refGenomeResourceMap) {
        this.refGenomeResourceMap = refGenomeResourceMap;
    }
}
