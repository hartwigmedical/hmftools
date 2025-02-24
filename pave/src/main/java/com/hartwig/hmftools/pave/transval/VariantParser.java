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

        ProteinVariant buildStopGained()
        {
            return new StopGained(geneData, transcriptData, aminoAcidsSequence, refRange);
        }

        ProteinVariant buildStartLost()
        {
            return new StartLost(geneData, transcriptData, aminoAcidsSequence, refRange);
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
        String[] geneVar = extractGeneAndVariant(expression);
        String variantDescription = geneVar[0];
        if(geneVar[1].startsWith("p."))
        {
            variantDescription = geneVar[1].substring(2);
        }
        return parseGeneVariant(geneVar[0], variantDescription);
    }

    public ProteinVariant parseGeneVariant(@NotNull String gene, @NotNull String variant)
    {
        if(variant.contains("delins"))
        {
            String[] refAltParts = variant.split("delins");
            return parseFactory(gene, variant).buildDeletionInsertion(refAltParts[1]);
        }
        if(variant.endsWith("del"))
        {
            return parseFactory(gene, variant).buildDeletion();
        }
        if(variant.endsWith("dup"))
        {
            return parseFactory(gene, variant).buildDuplication();
        }
        if(variant.contains("ins"))
        {
            String[] refAltParts = variant.split("ins");
            return parseFactory(gene, variant).buildInsertion(refAltParts[1]);
        }
        if(variant.endsWith("fs"))
        {
            return parseFactory(gene, variant).buildFrameshift();
        }
        if(variant.endsWith("*"))
        {
            return parseFactory(gene, variant).buildStopGained();
        }
        if(variant.endsWith("?"))
        {
            return parseFactory(gene, variant).buildStartLost();
        }
        GeneData geneData = lookupGene(gene);
        return buildSingleAminoAcidVariant(variant, geneData);
    }

    private VariantFactory parseFactory(@NotNull String gene, @NotNull String variantDescription)
    {
        GeneData geneData = lookupGene(gene);
        String rangeDescription = variantDescription.substring(0, variantDescription.length() - 3);
        AminoAcidRange refRange = getAminoAcidRange(variantDescription, rangeDescription);
        TranscriptData transcriptData = getApplicableTranscript(geneData, refRange);
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
    private SingleAminoAcidVariant buildSingleAminoAcidVariant(@NotNull final String variantDescription, final GeneData geneData)
    {
        Pattern variationPattern = Pattern.compile("^" + BLANK_OR_AA_GROUP + NAT + SINGLE_AA_GROUP + "$");
        final Matcher matcher = Checks.matchPattern(variationPattern, variantDescription);

        int position = Integer.parseInt(matcher.group(3));
        AminoAcidSpecification refSpec = new AminoAcidSpecification(position, matcher.group(1));
        AminoAcidSpecification altSpec = new AminoAcidSpecification(position, matcher.group(4));
        TranscriptData transcriptData = getApplicableTranscript(geneData, refSpec);
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

    public TranscriptData getApplicableTranscript(final GeneData geneData, final TranscriptFilter refFilter)
    {
        List<TranscriptData> allTranscripts = mEnsemblCache.getTranscripts(geneData.GeneId);
        List<TranscriptData> matchingRef = filter(allTranscripts, refFilter);
        if(matchingRef.isEmpty())
        {
            String msg = format("No transcript found for gene: %s matching ref: %s", geneData.GeneId, refFilter);
            throw new IllegalStateException(msg);
        }
        if(matchingRef.size() > 1)
        {
            return matchingRef.stream()
                    .filter(transcriptData -> transcriptData.IsCanonical)
                    .findFirst()
                    .orElseThrow(() ->
                    {
                        String msg =
                                format("No canonical transcript, but multiple non-canonical transcripts, found for gene: %s, ref: %s", geneData.GeneId, refFilter);
                        return new IllegalStateException(msg);
                    });
        }
        return matchingRef.get(0);
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
