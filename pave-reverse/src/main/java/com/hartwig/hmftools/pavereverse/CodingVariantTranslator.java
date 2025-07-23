package com.hartwig.hmftools.pavereverse;

import java.util.HashMap;
import java.util.List;
import java.util.Map;

import com.hartwig.hmftools.common.ensemblcache.EnsemblDataCache;
import com.hartwig.hmftools.common.ensemblcache.EnsemblDataLoader;
import com.hartwig.hmftools.common.gene.TranscriptAminoAcids;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeInterface;

/**
 * Translates an HGVS coding variant
 */
public class CodingVariantTranslator
{
    public final RefGenomeInterface RefGenome;
    public final EnsemblDataCache EnsemblCache;
    public final Map<String, TranscriptAminoAcids> TranscriptAminoAcidsMap = new HashMap<>();

    public CodingVariantTranslator(EnsemblDataCache ensemblCache, String ensemblDataPath, RefGenomeInterface refGenome)
    {
        EnsemblCache = ensemblCache;
        RefGenome = refGenome;
        EnsemblCache.setRequiredData(true, true, true, false);
        EnsemblCache.load(false);
        EnsemblCache.createTranscriptIdMap();
        EnsemblDataLoader.loadTranscriptAminoAcidData(ensemblDataPath, TranscriptAminoAcidsMap, List.of(), false);
    }

    public BaseSequenceChange translateHgvsCodingVariant(String gene, String transcript, String variant)
    {
        return null;
    }
}
