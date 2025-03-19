package com.hartwig.hmftools.pavereverse.parse;

import java.util.Map;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import com.google.common.annotations.VisibleForTesting;
import com.hartwig.hmftools.common.ensemblcache.EnsemblDataCache;
import com.hartwig.hmftools.common.gene.GeneData;
import com.hartwig.hmftools.common.gene.TranscriptAminoAcids;
import com.hartwig.hmftools.common.gene.TranscriptData;
import com.hartwig.hmftools.pavereverse.dna.DnaVariant;
import com.hartwig.hmftools.pavereverse.dna.HgvsAddress;
import com.hartwig.hmftools.pavereverse.dna.InExon;
import com.hartwig.hmftools.pavereverse.dna.InExonDownstreamOfCodingEnd;
import com.hartwig.hmftools.pavereverse.dna.InExonUpstreamOfCodingStart;
import com.hartwig.hmftools.pavereverse.dna.InIntronAfterExon;
import com.hartwig.hmftools.pavereverse.dna.InIntronBeforeExon;
import com.hartwig.hmftools.pavereverse.dna.InIntronDownstreamOfCodingEnd;
import com.hartwig.hmftools.pavereverse.dna.InIntronUpstreamOfCodingStart;

import org.apache.commons.lang3.NotImplementedException;

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
        String trimmedVariant = trimInitialCdot(variant);
        if(trimmedVariant.contains(">"))
        {
            Pattern roughStructurePattern = Pattern.compile("([+\\-*\\d]+)([ACGT]+)>([ACGT]+)");
            Matcher roughStructureMatcher = roughStructurePattern.matcher(trimmedVariant);
            boolean roughStructureFound = roughStructureMatcher.matches();

            if(!roughStructureFound)
            {
                throw new IllegalArgumentException("Variant not parsed: " + variant);
            }
            String positionInfo = roughStructureMatcher.group(1);
            HgvsAddress address = parseAddress(positionInfo);
            String ref = roughStructureMatcher.group(2);
            String alt = roughStructureMatcher.group(3);
            return new DnaVariant(geneData, transcriptData, address, ref, alt);
        }
        return null;
    }

    @VisibleForTesting
    static HgvsAddress parseAddress(String address)
    {
        // 67, -1, *1, 87+1, 88-1, -85+1, -84-3, *38-3, *23+2
        Pattern p = Pattern.compile("([\\-*]?)(\\d*?)([\\-+]?)(\\d+)");
        Matcher matcher = p.matcher(address);
        boolean matches = matcher.matches();
        if(!matches)
        {
            throw new IllegalArgumentException("Address not parsed: " + address);
        }
        String utrMarker = groupOrBlank(matcher, 1);
        String exonicBasePositionMarker = groupOrBlank(matcher, 2);
        String beforeAfterMarker = groupOrBlank(matcher, 3);
        int position = Integer.parseInt(groupOrBlank(matcher, 4));
        if("-".equals(utrMarker))
        {
            if(exonicBasePositionMarker.isBlank())
            {
                return new InExonUpstreamOfCodingStart(-1 * position);
            }
            else
            {
                int exonicBasePosition = -1 * Integer.parseInt(exonicBasePositionMarker);
                int relativePosition = plusOrMinus(beforeAfterMarker) * position;
                return new InIntronUpstreamOfCodingStart(exonicBasePosition, relativePosition);
            }
        }
        else if("*".equals(utrMarker))
        {
            if(exonicBasePositionMarker.isBlank())
            {
                return new InExonDownstreamOfCodingEnd(position);
            }
            else
            {
                int exonicBasePosition = Integer.parseInt(exonicBasePositionMarker);
                int relativePosition = plusOrMinus(beforeAfterMarker) * position;
                return new InIntronDownstreamOfCodingEnd(exonicBasePosition, relativePosition);
            }
        }
        else
        {
            if(exonicBasePositionMarker.isBlank())
            {
                return new InExon(position);
            }
            else
            {
                int exonBasePosition = Integer.parseInt(exonicBasePositionMarker);
                if(beforeAfterMarker.equals("+"))
                {
                    return new InIntronAfterExon(exonBasePosition, position);
                }
                else if(beforeAfterMarker.equals("-"))
                {
                    return new InIntronBeforeExon(exonBasePosition, position);
                }
                throw new NotImplementedException(address);
            }
        }
    }

    private static String groupOrBlank(Matcher matcher, int group)
    {
        String result = matcher.group(group);
        return result == null ? "" : result;
    }

    private static String trimInitialCdot(String variant)
    {
        if(variant.startsWith("c."))
        {
            return variant.substring(2);
        }
        return variant;
    }

    private static int plusOrMinus(String s)
    {
        if(s.equals("+"))
        {
            return 1;
        }
        else if(s.equals("-"))
        {
            return -1;
        }
        else
        {
            throw new IllegalArgumentException("Should be + or -:" + s);
        }
    }
}
