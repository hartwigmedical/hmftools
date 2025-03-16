package com.hartwig.hmftools.pavereverse.parse;

import static com.hartwig.hmftools.pavereverse.util.Checks.HGVS_FORMAT_REQUIRED;

import java.util.Map;

import com.hartwig.hmftools.common.ensemblcache.EnsemblDataCache;
import com.hartwig.hmftools.common.gene.GeneData;
import com.hartwig.hmftools.common.gene.TranscriptAminoAcids;

public class VariantParser
{
    final EnsemblDataCache EnsemblCache;
    final Map<String, TranscriptAminoAcids> TranscriptAminoAcidsMap;
    final String NAT = "(\\d+)";

    public VariantParser(final EnsemblDataCache ensemblDataCache, final Map<String, TranscriptAminoAcids> transcriptAminoAcidsMap)
    {
        EnsemblCache = ensemblDataCache;
        TranscriptAminoAcidsMap = transcriptAminoAcidsMap;
    }

    GeneData lookupGene(final String reference)
    {
        GeneData geneData = EnsemblCache.getGeneDataByName(reference);
        if(geneData == null)
        {
            geneData = EnsemblCache.getGeneDataById(reference);
            if(geneData == null)
            {
                throw new IllegalArgumentException(reference + " is not a known gene");
            }
        }
        return geneData;
    }

    static String[] extractGeneAndVariant(final String input)
    {
        String[] geneVar = input.split(":");
        if(geneVar.length != 2)
        {
            throw new IllegalArgumentException(HGVS_FORMAT_REQUIRED);
        }
        return geneVar;
    }
}
