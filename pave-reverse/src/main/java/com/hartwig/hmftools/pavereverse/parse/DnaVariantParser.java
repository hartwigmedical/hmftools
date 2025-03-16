package com.hartwig.hmftools.pavereverse.parse;

import java.util.Map;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import com.hartwig.hmftools.common.ensemblcache.EnsemblDataCache;
import com.hartwig.hmftools.common.gene.GeneData;
import com.hartwig.hmftools.common.gene.TranscriptAminoAcids;
import com.hartwig.hmftools.common.gene.TranscriptData;
import com.hartwig.hmftools.pavereverse.dna.DnaVariant;

public class DnaVariantParser extends VariantParser
{
    public DnaVariantParser(EnsemblDataCache ensemblDataCache, Map<String, TranscriptAminoAcids> transcriptAminoAcidsMap)
    {
        super(ensemblDataCache, transcriptAminoAcidsMap);
    }

    public DnaVariant parse(String gene, String transcriptId, String variant)
    {
        GeneData geneData = lookupGene(gene);
        TranscriptData transcriptData = EnsemblCache.getTranscriptData(geneData.GeneId, transcriptId);
        if(variant.contains(">"))
        {
            Pattern p = Pattern.compile("(\\d+)([ACGT]+)>([ACGT]+)");
            Matcher matcher = p.matcher(variant);
            boolean matches = matcher.find();
            if(!matches)
            {
                throw new IllegalArgumentException("Variant " + variant);
            }
            int position = Integer.parseInt(matcher.group(1));
            String ref = matcher.group(2);
            String alt = matcher.group(3);
            return new DnaVariant(geneData, transcriptData, position, ref, alt);
        }
        return null;
    }
}
