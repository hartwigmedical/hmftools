package com.hartwig.hmftools.pavereverse;

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
import com.hartwig.hmftools.pavereverse.variants.ProteinVariant;

/**
 * Calculates the possible base sequence changes that give rise to a protein-level change.
 */
public class ReversePave
{
    public final RefGenomeInterface RefGenome;
    public final EnsemblDataCache EnsemblCache;
    public final Map<String, TranscriptAminoAcids> TranscriptAminoAcidsMap = new HashMap<>();

    public ReversePave(EnsemblDataCache ensemblCache, String ensemblDataPath, RefGenomeInterface refGenome)
    {
        EnsemblCache = ensemblCache;
        RefGenome = refGenome;
        EnsemblCache.setRequiredData(true, true, true, false);
        EnsemblCache.load(false);
        EnsemblCache.createTranscriptIdMap();
        EnsemblDataLoader.loadTranscriptAminoAcidData(ensemblDataPath, TranscriptAminoAcidsMap, List.of(), false);
    }

    public ReversePave(File ensemblDataDir, RefGenomeVersion genomeVersion, RefGenomeInterface refGenome)
    {
        this(new EnsemblDataCache(ensemblDataDir.getAbsolutePath(), genomeVersion), ensemblDataDir.getAbsolutePath(), refGenome);
    }

    public VariantParser variationParser()
    {
        return new VariantParser(EnsemblCache, TranscriptAminoAcidsMap);
    }

    public BaseSequenceVariants calculateVariant(String proteinVariant)
    {
        ProteinVariant variant = variationParser().parse(proteinVariant);
        return variant.calculateVariant(RefGenome);
    }

    public BaseSequenceVariants calculateVariant(String gene, String proteinVariant)
    {
        ProteinVariant variant = variationParser().parseGeneVariant(gene, proteinVariant);
        return variant.calculateVariant(RefGenome);
    }

    public BaseSequenceVariants calculateVariant(String gene, String transcriptId, String proteinVariant)
    {
        ProteinVariant variant = variationParser().parseGeneVariant(gene, transcriptId, proteinVariant);
        return variant.calculateVariant(RefGenome);
    }

    public BaseSequenceVariants calculateVariantAllowMultipleNonCanonicalTranscriptMatches(String gene, String proteinVariant)
    {
        Set<ProteinVariant> allMatchingVariants = variationParser().parseGeneVariants(gene, proteinVariant);
        Map<String, BaseSequenceVariants> transcriptIdToVariant = new HashMap<>();
        allMatchingVariants.forEach(variant ->
        {
            try
            {
                transcriptIdToVariant.put(variant.Transcript.TransName, variant.calculateVariant(RefGenome));
            }
            catch(IllegalArgumentException e)
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
            transcriptIdToVariant.values().forEach(variant ->
            {
                if(!variant.changes().equals(example.changes()))
                {
                    String name1 = variant.Transcript.TransName;
                    String name2 = example.Transcript.TransName;
                    String msg = format("Transcripts %s and %s both match %s but produce different hotspots", name1, name2, proteinVariant);
                    throw new IllegalArgumentException(msg);
                }
            });
            return example;
        }
    }
}
