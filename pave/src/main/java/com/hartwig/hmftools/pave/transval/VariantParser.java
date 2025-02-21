package com.hartwig.hmftools.pave.transval;

import static java.lang.String.format;

import static com.hartwig.hmftools.pave.transval.Checks.HGVS_FORMAT_REQUIRED;

import java.util.List;
import java.util.Map;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import java.util.stream.Collectors;

import com.hartwig.hmftools.common.ensemblcache.EnsemblDataCache;
import com.hartwig.hmftools.common.gene.GeneData;
import com.hartwig.hmftools.common.gene.TranscriptAminoAcids;
import com.hartwig.hmftools.common.gene.TranscriptData;

import org.jetbrains.annotations.NotNull;

class VariantParser
{
    private final EnsemblDataCache mEnsemblCache;
    private final Map<String, TranscriptAminoAcids> mTranscriptAminoAcidsMap;
    final String SINGLE_AA = "[A-Z]([a-z][a-z])?+";
    final String SINGLE_AA_GROUP = "(" + SINGLE_AA + ")";
    final String BLANK_OR_AA_GROUP = "(" + SINGLE_AA + ")?+";
    final String NAT = "(\\d+)";

    private static class VariantFactory
    {
        @NotNull
        private final GeneData geneData;
        @NotNull
        private final AminoAcidRange refRange;
        @NotNull
        private final TranscriptData transcriptData;
        @NotNull
        private final TranscriptAminoAcids aminoAcidsSequence;


        private VariantFactory(@NotNull final GeneData geneData, @NotNull final AminoAcidRange refRange,
                @NotNull final TranscriptData transcriptData, @NotNull final TranscriptAminoAcids aminoAcidsSequence)
        {
            this.geneData = geneData;
            this.refRange = refRange;
            this.transcriptData = transcriptData;
            this.aminoAcidsSequence = aminoAcidsSequence;
        }

        ProteinVariant buildDuplication()
        {
            return new Duplication(geneData, transcriptData, aminoAcidsSequence, refRange);
        }

        ProteinVariant buildFrameshift()
        {
            return new Frameshift(geneData, transcriptData, aminoAcidsSequence, refRange);
        }

        ProteinVariant buildDeletion()
        {
            return new Deletion(geneData, transcriptData, aminoAcidsSequence, refRange);
        }

        ProteinVariant buildDeletionInsertion(String altSequence)
        {
            AminoAcidSequence altAminoAcids = AminoAcidSequence.parse(altSequence);
            return new DeletionInsertion(geneData, transcriptData, aminoAcidsSequence, refRange, altAminoAcids);
        }

        ProteinVariant buildInsertion(String inserted)
        {
            AminoAcidSequence altAminoAcids = AminoAcidSequence.parse(inserted);
            return new Insertion(geneData, transcriptData, aminoAcidsSequence, refRange, altAminoAcids);
        }
    }

    public VariantParser(final EnsemblDataCache ensemblDataCache, final Map<String, TranscriptAminoAcids> transcriptAminoAcidsMap)
    {
        mEnsemblCache = ensemblDataCache;
        mTranscriptAminoAcidsMap = transcriptAminoAcidsMap;
    }

    public ProteinVariant parse(@NotNull String expression)
    {
        if(expression.contains("delins"))
        {
            return parseDeletionInsertion(expression);
        }
        if(expression.contains("del"))
        {
            return parseDeletion(expression);
        }
        return parseSingleAminoAcidVariant(expression);
    }

    public ProteinVariant parseVariantForGene(@NotNull String gene, @NotNull String variant)
    {
        if(variant.contains("delins"))
        {
            return parseDeletionInsertion(gene, variant);
        }
        if(variant.endsWith("del"))
        {
            return parseDeletion(gene, variant);
        }
        if(variant.endsWith("dup"))
        {
            return parseDuplication(gene, variant);
        }
        if(variant.contains("ins"))
        {
            return parseInsertion(gene, variant);
        }
        if(variant.endsWith("fs"))
        {
            return parseFactory(gene, variant).buildFrameshift();
        }
        return parseSingleAminoAcidVariant(gene, variant);
    }

    public SingleAminoAcidVariant parseSingleAminoAcidVariant(@NotNull String gene, @NotNull String variantDescription)
    {
        GeneData geneData = lookupGene(gene);
        return buildSingleAminoAcidVariant(variantDescription, geneData);
    }

    public SingleAminoAcidVariant parseSingleAminoAcidVariant(String input)
    {
        String[] geneVar = extractGeneAndVariant(input);
        String variantDescription = extractDescription(geneVar)[1];
        GeneData geneData = lookupGene(geneVar[0]);
        return buildSingleAminoAcidVariant(variantDescription, geneData);
    }

    public ProteinVariant parseDeletionInsertion(@NotNull String gene, @NotNull String variantDescription)
    {
        String[] refAltParts = variantDescription.split("delins");
        return parseFactory(gene, variantDescription).buildDeletionInsertion(refAltParts[1]);
    }

    public ProteinVariant parseDeletion(@NotNull String gene, @NotNull String variantDescription)
    {
        return parseFactory(gene, variantDescription).buildDeletion();
    }

    public ProteinVariant parseDeletionInsertion(String input)
    {
        String[] geneVar = extractGeneAndVariant(input);
        return parseDeletionInsertion(geneVar[0], geneVar[1]);
    }

    public ProteinVariant parseDeletion(String input)
    {
        String[] geneVar = extractGeneAndVariant(input);
        return parseFactory(geneVar[0], geneVar[1]).buildDeletion();
    }

    public ProteinVariant parseDuplication(@NotNull String gene, @NotNull String variantDescription)
    {
        return parseFactory(gene, variantDescription).buildDuplication();
    }

    public ProteinVariant parseInsertion(@NotNull String gene, @NotNull String description)
    {
        String[] refAltParts = description.split("ins");
        return parseFactory(gene, description).buildInsertion(refAltParts[1]);
    }

    private VariantFactory parseFactory(@NotNull String gene, @NotNull String variantDescription)
    {
        GeneData geneData = lookupGene(gene);
        String rangeDescription = variantDescription.substring(0, variantDescription.length() - 3);
        AminoAcidRange refRange = getAminoAcidRange(variantDescription, rangeDescription);
        TranscriptData transcriptData = getApplicableTranscript(geneData, refRange, new PassThroughFilter());
        TranscriptAminoAcids aminoAcidsSequence = lookupTranscriptAminoAcids(transcriptData, false);
        return new VariantFactory(geneData, refRange, transcriptData, aminoAcidsSequence);
    }

    @NotNull
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

    @NotNull
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

    @NotNull
    private SingleAminoAcidVariant buildSingleAminoAcidVariant(@NotNull final String variantDescription,
            final GeneData geneData)
    {
        Pattern variationPattern = Pattern.compile("^" + BLANK_OR_AA_GROUP + NAT + SINGLE_AA_GROUP + "$");
        final Matcher matcher = Checks.matchPattern(variationPattern, variantDescription);

        int position = Integer.parseInt(matcher.group(3));
        AminoAcidSpecification refSpec = new AminoAcidSpecification(position, matcher.group(1));
        AminoAcidSpecification altSpec = new AminoAcidSpecification(position, matcher.group(4));
        TranscriptData transcriptData = getApplicableTranscript(geneData, refSpec, new Negation(altSpec));
        TranscriptAminoAcids aminoAcidsSequence = lookupTranscriptAminoAcids(transcriptData, false);
        return new SingleAminoAcidVariant(geneData, transcriptData, aminoAcidsSequence, altSpec);
    }

    private TranscriptAminoAcids lookupTranscriptAminoAcids(final TranscriptData transcriptData, final boolean allowNullReturn)
    {
        TranscriptAminoAcids aminoAcidsSequence = mTranscriptAminoAcidsMap.get(transcriptData.TransName);
        if(aminoAcidsSequence == null && !allowNullReturn)
        {
            throw new IllegalStateException("No amino acid sequence found for " + transcriptData.TransName);
        }
        return aminoAcidsSequence;
    }

    private List<TranscriptData> filter(List<TranscriptData> transcriptDataList, TranscriptFilter filter)
    {
        return transcriptDataList.stream().filter(transcriptData ->
        {
            TranscriptAminoAcids aminoAcidsSequence = lookupTranscriptAminoAcids(transcriptData, true);
            return aminoAcidsSequence != null && filter.applies(aminoAcidsSequence);
        }).collect(Collectors.toList());
    }

    public TranscriptData getApplicableTranscript(final GeneData geneData,
            final TranscriptFilter refFilter,
            final TranscriptFilter altFilter)
    {
        List<TranscriptData> allTranscripts = mEnsemblCache.getTranscripts(geneData.GeneId);
        List<TranscriptData> matchingRef = filter(allTranscripts, refFilter);
        if(matchingRef.isEmpty())
        {
            String msg = format("No transcript found for gene: %s matching ref: %s", geneData.GeneId, refFilter);
            throw new IllegalStateException(msg);
        }
        List<TranscriptData> possibilities = filter(matchingRef, altFilter);
        if(possibilities.isEmpty())
        {
            String msg =
                    format("No transcript found matching ref but not alt for gene: %s, ref: %s, alt: %s", geneData.GeneId, refFilter, altFilter);
            throw new IllegalStateException(msg);
        }
        if(possibilities.size() > 1)
        {
            return possibilities.stream()
                    .filter(transcriptData -> transcriptData.IsCanonical)
                    .findFirst()
                    .orElseThrow(() ->
                    {
                        String msg =
                                format("No canonical transcript, but multiple non-canonical transcripts, found for gene: %s, ref: %s, alt: %s", geneData.GeneId, refFilter, altFilter);
                        return new IllegalStateException(msg);
                    });
        }
        return possibilities.get(0);
    }

    @NotNull
    private GeneData lookupGene(final String reference)
    {
        GeneData geneData = mEnsemblCache.getGeneDataByName(reference);
        if(geneData == null)
        {
            geneData = mEnsemblCache.getGeneDataById(reference);
            if(geneData == null)
            {
                throw new IllegalArgumentException(reference + " is not a known gene");
            }
        }
        return geneData;
    }

    private static String[] extractDescription(final String[] geneVar)
    {
        String[] descriptionParts = geneVar[1].split("\\.");
        if(descriptionParts.length != 2)
        {
            throw new IllegalArgumentException(HGVS_FORMAT_REQUIRED);
        }
        if(!descriptionParts[0].equals("p"))
        {
            throw new IllegalArgumentException(HGVS_FORMAT_REQUIRED);
        }
        return descriptionParts;
    }

    private static String[] extractGeneAndVariant(final String input)
    {
        String[] geneVar = input.split(":");
        if(geneVar.length != 2)
        {
            throw new IllegalArgumentException(HGVS_FORMAT_REQUIRED);
        }
        return geneVar;
    }
}
