package com.hartwig.hmftools.pave.reverse;

import static java.lang.String.format;

import java.io.File;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;

import com.hartwig.hmftools.common.ensemblcache.EnsemblDataCache;
import com.hartwig.hmftools.common.ensemblcache.EnsemblDataLoader;
import com.hartwig.hmftools.common.gene.TranscriptAminoAcids;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeInterface;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion;

/**
 * Calculates the possible base sequence changes that give rise to a protein-level change.
 */
public class ReversePave
{
    public final RefGenomeInterface mRefGenome;
    public final EnsemblDataCache mEnsemblCache;
    final Map<String, TranscriptAminoAcids> mTranscriptAminoAcidsMap = new HashMap<>();

    public ReversePave(final EnsemblDataCache ensemblCache, final RefGenomeInterface refGenome)
    {
        mEnsemblCache = ensemblCache;
        mRefGenome = refGenome;
        mEnsemblCache.setRequiredData(true, true, true, false);
        mEnsemblCache.load(false);
        mEnsemblCache.createTranscriptIdMap();
        EnsemblDataLoader.loadTranscriptAminoAcidData(mEnsemblCache.dataPath(), mTranscriptAminoAcidsMap, List.of(), false);
    }

    public ReversePave(final File ensemblDataDir, final RefGenomeVersion genomeVersion, final RefGenomeInterface refGenome)
    {
        this(new EnsemblDataCache(ensemblDataDir.getAbsolutePath(), genomeVersion), refGenome);
    }

    VariantParser variationParser()
    {
        return new VariantParser(mEnsemblCache, mTranscriptAminoAcidsMap);
    }

    public BaseSequenceVariants calculateVariant(String proteinVariant)
    {
        ProteinVariant variant = variationParser().parse(proteinVariant);
        return variant.calculateVariant(mRefGenome);
    }

    public BaseSequenceVariants calculateVariant(String gene, String proteinVariant)
    {
        ProteinVariant variant = variationParser().parseGeneVariant(gene, proteinVariant);
        return variant.calculateVariant(mRefGenome);
    }

    public BaseSequenceVariants calculateVariant(String gene, String transcriptId, String proteinVariant)
    {
        ProteinVariant variant = variationParser().parseGeneVariant(gene, transcriptId, proteinVariant);
        return variant.calculateVariant(mRefGenome);
    }

    public BaseSequenceVariants calculateVariantAllowMultipleNonCanonicalTranscriptMatches(String gene, String proteinVariant)
    {
        Set<ProteinVariant> allMatchingVariants = variationParser().parseGeneVariants(gene, proteinVariant);
        Map<String, BaseSequenceVariants> transcriptIdToVariant  = new HashMap<>();
        allMatchingVariants.forEach(variant -> {
            try
            {
                transcriptIdToVariant.put(variant.mTranscript.TransName, variant.calculateVariant(mRefGenome));
            }
            catch (IllegalArgumentException e)
            {
                // eg incomplete transcript... just ignore it.
            }
        });
        if(transcriptIdToVariant.isEmpty())
        {
            // Maybe the only transcripts had problems.
            throw new IllegalArgumentException("No variant could be calculated for " + gene + " " + proteinVariant);
        }
        else if(transcriptIdToVariant.size() == 1)
        {
            return transcriptIdToVariant.values().iterator().next();
        }
        else
        {
            BaseSequenceVariants example = transcriptIdToVariant.values().iterator().next();
            transcriptIdToVariant.values().forEach(variant -> {
               if(!variant.changes().equals(example.changes()))
               {
                   String name1 = variant.mTranscript.TransName;
                   String name2 = example.mTranscript.TransName;
                   String msg = format("Transcripts %s and %s both match %s but produce different hotspots", name1, name2, proteinVariant);
                   throw new IllegalArgumentException(msg);
               }
            });
            return example;
        }
    }
}
