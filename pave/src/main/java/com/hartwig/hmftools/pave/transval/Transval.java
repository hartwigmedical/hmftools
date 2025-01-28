package com.hartwig.hmftools.pave.transval;

import java.io.File;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import com.hartwig.hmftools.common.ensemblcache.EnsemblDataCache;
import com.hartwig.hmftools.common.ensemblcache.EnsemblDataLoader;
import com.hartwig.hmftools.common.gene.TranscriptAminoAcids;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeInterface;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion;

public class Transval
{
    public final RefGenomeInterface mRefGenome;
    public final EnsemblDataCache mEnsemblCache;
    private final Map<String, TranscriptAminoAcids> mTranscriptAminoAcidsMap = new HashMap<>();

    public Transval(final File ensemblDataDir, final RefGenomeInterface refGenomeVersion)
    {
        this.mEnsemblCache = new EnsemblDataCache(ensemblDataDir.getAbsolutePath(), RefGenomeVersion.V38);
        this.mRefGenome = refGenomeVersion;
        mEnsemblCache.setRequiredData(true, true, true, false);
        mEnsemblCache.load(false);
        mEnsemblCache.createTranscriptIdMap();
        EnsemblDataLoader.loadTranscriptAminoAcidData(ensemblDataDir, mTranscriptAminoAcidsMap, List.of(), false);
    }

    public VariationParser variationParser()
    {
        return new VariationParser(mEnsemblCache, mTranscriptAminoAcidsMap);
    }

    public TransvalVariant calculateVariant(String proteinVariant)
    {
        ProteinVariant variant = variationParser().parse(proteinVariant);
        return variant.calculateVariant(mRefGenome);
    }
}
