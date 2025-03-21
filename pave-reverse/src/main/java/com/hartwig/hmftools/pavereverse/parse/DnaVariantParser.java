package com.hartwig.hmftools.pavereverse.parse;

import java.util.Map;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import com.google.common.base.Preconditions;
import com.hartwig.hmftools.common.ensemblcache.EnsemblDataCache;
import com.hartwig.hmftools.common.gene.GeneData;
import com.hartwig.hmftools.common.gene.TranscriptAminoAcids;
import com.hartwig.hmftools.common.gene.TranscriptData;
import com.hartwig.hmftools.pavereverse.dna.DeletionVariant;
import com.hartwig.hmftools.pavereverse.dna.DnaVariant;
import com.hartwig.hmftools.pavereverse.dna.DuplicationVariant;
import com.hartwig.hmftools.pavereverse.dna.HgvsAddress;
import com.hartwig.hmftools.pavereverse.dna.InExon;
import com.hartwig.hmftools.pavereverse.dna.InExonDownstreamOfCodingEnd;
import com.hartwig.hmftools.pavereverse.dna.InExonUpstreamOfCodingStart;
import com.hartwig.hmftools.pavereverse.dna.InIntronAfterExon;
import com.hartwig.hmftools.pavereverse.dna.InIntronBeforeExon;
import com.hartwig.hmftools.pavereverse.dna.InIntronDownstreamOfCodingEnd;
import com.hartwig.hmftools.pavereverse.dna.InIntronUpstreamOfCodingStart;
import com.hartwig.hmftools.pavereverse.dna.SubstitutionVariant;

import org.apache.commons.lang3.NotImplementedException;
import org.apache.commons.lang3.tuple.ImmutablePair;
import org.apache.commons.lang3.tuple.Pair;

public class DnaVariantParser extends VariantParser
{
    private static final String ADDRESS_RANGE_CHARS = "([+\\-*\\d_]+)";
    private static final String OPTIONAL_BASES = "([ACGT]*)";

    private interface VariantBuilder
    {
        DnaVariant build(GeneData geneData,
                TranscriptData transcriptData,
                Pair<HgvsAddress, HgvsAddress> addresses,
                String ref,
                String alt);
    }

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
            return getSubstitutionVariant(trimmedVariant, geneData, transcriptData);
        }
        else if(trimmedVariant.contains("dup"))
        {
            return getDupVariant(trimmedVariant, geneData, transcriptData);
        }
        else if(trimmedVariant.contains("del"))
        {
            return getDelVariant(trimmedVariant, geneData, transcriptData);
        }
        return null;
    }

    private static DnaVariant buildVariant(String variant,
            String variantSpecificPattern,
            GeneData geneData,
            TranscriptData transcriptData,
            VariantBuilder variantBuilder)
    {
        Pattern roughStructurePattern = Pattern.compile(ADDRESS_RANGE_CHARS + variantSpecificPattern);
        Matcher roughStructureMatcher = roughStructurePattern.matcher(variant);
        boolean roughStructureFound = roughStructureMatcher.matches();
        if(!roughStructureFound)
        {
            throw new IllegalArgumentException("Variant not parsed: " + variant);
        }
        String positionInfo = roughStructureMatcher.group(1);
        Pair<HgvsAddress, HgvsAddress> addresses = parseAddresses(positionInfo);
        String ref = groupOrBlank(roughStructureMatcher, 2);
        String alt = groupOrBlank(roughStructureMatcher, 3);
        return variantBuilder.build(geneData, transcriptData, addresses, ref, alt);
    }

    private static DnaVariant getDelVariant(String variant, GeneData geneData, TranscriptData transcriptData)
    {
        return buildVariant(variant,
                "del" + OPTIONAL_BASES,
                geneData,
                transcriptData,
                (geneData1, transcriptData1, addresses, ref, alt) -> new DeletionVariant(geneData1, transcriptData1, addresses.getLeft(), addresses.getRight(), ref));
    }

    private static DnaVariant getDupVariant(String variant, GeneData geneData, TranscriptData transcriptData)
    {
        return buildVariant(variant,
                "dup" + OPTIONAL_BASES,
                geneData,
                transcriptData,
                (geneData1, transcriptData1, addresses, ref, alt) -> new DuplicationVariant(geneData1, transcriptData1, addresses.getLeft(), addresses.getRight(), ref));
    }

    private static DnaVariant getSubstitutionVariant(String variant, GeneData geneData, TranscriptData transcriptData)
    {
        return buildVariant(variant,
                "([ACGT]+)>([ACGT]+)",
                geneData,
                transcriptData,
                (geneData1, transcriptData1, addresses, ref, alt) -> new SubstitutionVariant(geneData, transcriptData, addresses.getLeft(), ref, alt));
    }

    private static Pair<HgvsAddress, HgvsAddress> parseAddresses(String addresses)
    {
        if(addresses.contains("_"))
        {
            String[] parts = addresses.split("_");
            Preconditions.checkArgument(parts.length == 2, "Invalid address: " + addresses);
            return new ImmutablePair<>(parseAddress(parts[0]), parseAddress(parts[1]));
        }
        HgvsAddress address = parseAddress(addresses);
        return new ImmutablePair<>(address, address);
    }

    private static HgvsAddress parseAddress(String address)
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
        if(group > matcher.groupCount())
        {
            return "";
        }
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
