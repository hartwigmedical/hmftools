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
    public final RefGenomeInterface mRefGenomeVersion;
    public final EnsemblDataCache mEnsemblCache;
    private final Map<String, TranscriptAminoAcids> transcriptAminoAcidsMap = new HashMap<>();

    public Transval(final File ensemblDataDir, final RefGenomeInterface refGenomeVersion)
    {
        this.mEnsemblCache = new EnsemblDataCache(ensemblDataDir.getAbsolutePath(), RefGenomeVersion.V38);
        this.mRefGenomeVersion = refGenomeVersion;
        mEnsemblCache.setRequiredData(true, true, true, true);
        mEnsemblCache.load(false);
        mEnsemblCache.createTranscriptIdMap();
        EnsemblDataLoader.loadTranscriptAminoAcidData(ensemblDataDir, transcriptAminoAcidsMap, List.of(), false);
    }

    public VariationParser variationParser()
    {
        return new VariationParser(mEnsemblCache, transcriptAminoAcidsMap);
    }

    public TransvalSNV calculateSNV(String proteinVariant)
    {
        SingleAminoAcidVariant variant = variationParser().parse(proteinVariant);

        return new TransvalSNV(
                variant.Transcript.TransName,
                variant.Gene.Chromosome,
                variant.Position,
                !variant.codonIsInSingleExon(),
                variant.referenceAminoAcid(), // todo fix
                variant.variantAminoAcid(), // fix
                variant.referenceCodon(mRefGenomeVersion),
                List.of()
        );
    }
}
