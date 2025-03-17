package com.hartwig.hmftools.pavereverse.parse;

import java.util.Map;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import com.hartwig.hmftools.common.ensemblcache.EnsemblDataCache;
import com.hartwig.hmftools.common.gene.GeneData;
import com.hartwig.hmftools.common.gene.TranscriptAminoAcids;
import com.hartwig.hmftools.common.gene.TranscriptData;
import com.hartwig.hmftools.pavereverse.dna.DnaVariant;
import com.hartwig.hmftools.pavereverse.dna.DownstreamOfCodingEndAddress;
import com.hartwig.hmftools.pavereverse.dna.HgvsAddress;
import com.hartwig.hmftools.pavereverse.dna.InExonAddress;
import com.hartwig.hmftools.pavereverse.dna.UpstreamOfCodingStartAddress;

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
            Pattern p = Pattern.compile("([-*])?(\\d+)([ACGT]+)>([ACGT]+)");
            Matcher matcher = p.matcher(variant);
            boolean matches = matcher.find();
            if(!matches)
            {
                throw new IllegalArgumentException("Variant " + variant);
            }
            String upstreamDownstreamMarker = matcher.group(1);
            int position = Integer.parseInt(matcher.group(2));
            String ref = matcher.group(3);
            String alt = matcher.group(4);
            HgvsAddress address;
            if("-".equals(upstreamDownstreamMarker))
            {
                address = new UpstreamOfCodingStartAddress(position);
            }
            else if("*".equals(upstreamDownstreamMarker))
            {
                address = new DownstreamOfCodingEndAddress(position);
            }
            else
            {
                address = new InExonAddress(position);
            }
            return new DnaVariant(geneData, transcriptData, address, ref, alt);
        }
        return null;
    }
}
