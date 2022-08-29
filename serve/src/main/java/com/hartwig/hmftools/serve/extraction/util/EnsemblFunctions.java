package com.hartwig.hmftools.serve.extraction.util;

import com.hartwig.hmftools.common.ensemblcache.EnsemblDataCache;
import com.hartwig.hmftools.common.gene.GeneData;
import com.hartwig.hmftools.common.gene.TranscriptData;

import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public final class EnsemblFunctions {

    private EnsemblFunctions() {
    }

    @Nullable
    public static HmfTranscriptRegion findCanonicalTranscript(@NotNull EnsemblDataCache ensemblDataCache, @NotNull String gene) {
        GeneData geneData = ensemblDataCache.getGeneDataByName(gene);
        if (geneData == null) {
            return null;
        }

        TranscriptData transcriptData = ensemblDataCache.getCanonicalTranscriptData(geneData.GeneId);
        if (transcriptData == null) {
            return null;
        }

        return HmfTranscriptRegionUtils.fromTranscript(geneData, transcriptData);
    }
}
