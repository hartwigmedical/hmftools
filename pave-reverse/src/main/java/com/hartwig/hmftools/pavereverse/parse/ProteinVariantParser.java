package com.hartwig.hmftools.pavereverse.parse;

import static com.hartwig.hmftools.pavereverse.util.Checks.HGVS_FORMAT_REQUIRED;

import java.util.Map;
import java.util.Set;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import com.google.common.base.Preconditions;
import com.hartwig.hmftools.common.ensemblcache.EnsemblDataCache;
import com.hartwig.hmftools.common.gene.GeneData;
import com.hartwig.hmftools.common.gene.TranscriptAminoAcids;
import com.hartwig.hmftools.pavereverse.aa.AminoAcidRange;
import com.hartwig.hmftools.pavereverse.aa.AminoAcidSpecification;
import com.hartwig.hmftools.pavereverse.protein.ProteinVariant;
import com.hartwig.hmftools.pavereverse.util.Checks;

public class ProteinVariantParser extends VariantParser
{
    final String SINGLE_AA = "[A-Z]([a-z][a-z])?+";
    final String SINGLE_AA_GROUP = "(" + SINGLE_AA + ")";
    final String BLANK_OR_AA_GROUP = "(" + SINGLE_AA + ")?+";

    public ProteinVariantParser(final EnsemblDataCache ensemblDataCache, final Map<String, TranscriptAminoAcids> transcriptAminoAcidsMap)
    {
        super(ensemblDataCache, transcriptAminoAcidsMap);
    }

    public Set<ProteinVariant> parseGeneVariants(String gene, String variant)
    {
        return parseGeneVariant(gene, variant, new FilterTranscripts(this, true));
    }

    public ProteinVariant parseGeneVariant(String gene, String variant)
    {
        Set<ProteinVariant> variants = parseGeneVariant(gene, variant, new FilterTranscripts(this, false));
        Preconditions.checkArgument(variants.size() == 1);
        return variants.iterator().next();
    }

    public ProteinVariant parseGeneVariant(String gene, String specificTranscript, String variant)
    {
        return parseGeneVariant(gene, trimInitialPdot(variant), new GetSpecificTranscript(this, specificTranscript)).iterator().next();
    }

    public ProteinVariant parse(String expression)
    {
        String[] geneVar = extractGeneAndVariant(expression);
        String variantDescription = trimInitialPdot(geneVar[1]);
        return parseGeneVariant(geneVar[0], variantDescription);
    }

    private static String trimInitialPdot(String variant)
    {
        if(variant.startsWith("p."))
        {
            return variant.substring(2);
        }
        return variant;
    }

    private Set<ProteinVariant> parseGeneVariant(String gene, String variant, TranscriptRetriever transcriptRetriever)
    {
        if(variant.contains("delins"))
        {
            String[] refAltParts = variant.split("delins");
            return parseFactory(gene, variant, transcriptRetriever).buildDeletionInsertion(refAltParts[1]);
        }
        if(variant.endsWith("del"))
        {
            return parseFactory(gene, variant, transcriptRetriever).buildDeletion();
        }
        if(variant.endsWith("dup"))
        {
            return parseFactory(gene, variant, transcriptRetriever).buildDuplication();
        }
        if(variant.contains("ins"))
        {
            String[] refAltParts = variant.split("ins");
            return parseFactory(gene, variant, transcriptRetriever).buildInsertion(refAltParts[1]);
        }
        if(variant.endsWith("fs"))
        {
            return parseFactory(gene, variant, transcriptRetriever).buildFrameshift();
        }
        if(variant.endsWith("*"))
        {
            return parseFactory(gene, variant, transcriptRetriever).buildStopGained();
        }
        if(variant.endsWith("?"))
        {
            return parseFactory(gene, variant, transcriptRetriever).buildStartLost();
        }
        Pattern variationPattern = Pattern.compile("^" + BLANK_OR_AA_GROUP + NAT + SINGLE_AA_GROUP + "$");
        final Matcher matcher = Checks.matchPattern(variationPattern, variant);
        int position = Integer.parseInt(matcher.group(3));
        AminoAcidSpecification altSpec = new AminoAcidSpecification(position, matcher.group(4));
        return parseFactory(gene, variant, transcriptRetriever).buildSingleAminoAcidVariants(altSpec);
    }

    private ProteinVariantFactory parseFactory(String gene, String variantDescription,
            final TranscriptRetriever transcriptRetriever)
    {
        GeneData geneData = lookupGene(gene);
        String rangeDescription = variantDescription.substring(0, variantDescription.length() - 3);
        AminoAcidRange refRange = getAminoAcidRange(variantDescription, rangeDescription);
        Set<ProteinTranscript> transcriptData = transcriptRetriever.getApplicableTranscripts(geneData, refRange);
        return new ProteinVariantFactory(geneData, refRange, transcriptData);
    }

    private static AminoAcidRange getAminoAcidRange(final String description, final String rangeDescription)
    {
        AminoAcidRange refRange;
        if(rangeDescription.contains("_"))
        {
            refRange = parseAminoAcidRange(rangeDescription);
        }
        else
        {
            AminoAcidSpecification refStart = AminoAcidSpecification.parse(description);
            refRange = new AminoAcidRange(refStart, refStart);
        }
        return refRange;
    }

    private static AminoAcidRange parseAminoAcidRange(final String refAltParts)
    {
        String[] refParts = refAltParts.split("_");
        if(refParts.length != 2)
        {
            throw new IllegalArgumentException(HGVS_FORMAT_REQUIRED);
        }
        AminoAcidSpecification refStart = AminoAcidSpecification.parse(refParts[0]);
        AminoAcidSpecification refStop = AminoAcidSpecification.parse(refParts[1]);
        return new AminoAcidRange(refStart, refStop);
    }
}
