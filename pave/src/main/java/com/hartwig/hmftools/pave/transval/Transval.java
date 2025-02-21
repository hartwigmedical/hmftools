package com.hartwig.hmftools.pave.transval;

import java.io.File;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import com.google.common.base.Preconditions;
import com.hartwig.hmftools.common.ensemblcache.EnsemblDataCache;
import com.hartwig.hmftools.common.ensemblcache.EnsemblDataLoader;
import com.hartwig.hmftools.common.gene.TranscriptAminoAcids;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeInterface;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion;

public class Transval
{
    public final RefGenomeInterface mRefGenome;
    public final EnsemblDataCache mEnsemblCache;
    final Map<String, TranscriptAminoAcids> mTranscriptAminoAcidsMap = new HashMap<>();

    public Transval(final File ensemblDataDir, final RefGenomeInterface refGenomeVersion)
    {
        Preconditions.checkArgument(ensemblDataDir.isDirectory());
        Preconditions.checkArgument(ensemblDataDir.exists());
        this.mEnsemblCache = new EnsemblDataCache(ensemblDataDir.getAbsolutePath(), RefGenomeVersion.V38);
        this.mRefGenome = refGenomeVersion;
        mEnsemblCache.setRequiredData(true, true, true, false);
        mEnsemblCache.load(false);
        mEnsemblCache.createTranscriptIdMap();
        EnsemblDataLoader.loadTranscriptAminoAcidData(ensemblDataDir, mTranscriptAminoAcidsMap, List.of(), false);
    }

    VariantParser variationParser()
    {
        return new VariantParser(mEnsemblCache, mTranscriptAminoAcidsMap);
    }

    public TransvalVariant calculateVariant(String proteinVariant)
    {
        ProteinVariant variant = variationParser().parse(proteinVariant);
        return variant.calculateVariant(mRefGenome);
    }

    public TransvalVariant calculateVariant(String gene, String proteinVariant)
    {
        ProteinVariant variant = variationParser().parseVariantForGene(gene, proteinVariant);
        return variant.calculateVariant(mRefGenome);
    }
}
