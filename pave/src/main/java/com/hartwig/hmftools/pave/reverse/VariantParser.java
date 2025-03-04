package com.hartwig.hmftools.pave.reverse;

import static java.lang.String.format;

import static com.hartwig.hmftools.pave.reverse.Checks.HGVS_FORMAT_REQUIRED;

import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import java.util.stream.Collectors;

import com.google.common.base.Preconditions;
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

    private static class ProteinTranscript
    {
        @NotNull
        final TranscriptAminoAcids mAminoAcids;
        @NotNull
        final TranscriptData mTranscriptData;

        private ProteinTranscript(@NotNull final TranscriptAminoAcids aminoAcids, @NotNull final TranscriptData transcriptData)
        {
            this.mAminoAcids = aminoAcids;
            this.mTranscriptData = transcriptData;
        }
    }

    private interface TranscriptRetriever
    {
        Set<ProteinTranscript> getApplicableTranscripts(final GeneData geneData, final TranscriptFilter refFilter);
    }

    private class FilterTranscripts implements TranscriptRetriever
    {
        final boolean returnVariantForEachNonCanonicalTranscript;

        private FilterTranscripts(final boolean returnVariantForEachNonCanonicalTranscript)
        {
            this.returnVariantForEachNonCanonicalTranscript = returnVariantForEachNonCanonicalTranscript;
        }

        @Override
        public Set<ProteinTranscript> getApplicableTranscripts(final GeneData geneData, final TranscriptFilter refFilter)
        {
            List<TranscriptData> allTranscripts = mEnsemblCache.getTranscripts(geneData.GeneId);
            Set<ProteinTranscript> matchingRef = filter(allTranscripts, refFilter);
            if(matchingRef.isEmpty())
            {
                String msg = format("No transcript found for gene: %s matching ref: %s", geneData.GeneId, refFilter);
                throw new IllegalArgumentException(msg);
            }
            ProteinTranscript canonical = matchingRef.stream()
                    .filter(pt -> pt.mTranscriptData.IsCanonical)
                    .findFirst()
                    .orElse(null);
            if(canonical != null)
            {
                return Set.of(canonical);
            }
            if(matchingRef.size() > 1 && !returnVariantForEachNonCanonicalTranscript)
            {
                String msg =
                        format("No canonical transcript, but multiple non-canonical transcripts, found for gene: %s, ref: %s", geneData.GeneId, refFilter);
                throw new IllegalArgumentException(msg);

            }
            return new HashSet<>(matchingRef);
        }

        private Set<ProteinTranscript> filter(List<TranscriptData> transcriptDataList, TranscriptFilter filter)
        {
            Set<ProteinTranscript> result = new HashSet<>();
            transcriptDataList.forEach(transcriptData ->
            {
                TranscriptAminoAcids aminoAcidsSequence = mTranscriptAminoAcidsMap.get(transcriptData.TransName);
                if(aminoAcidsSequence != null && filter.applies(aminoAcidsSequence))
                {
                    result.add(new ProteinTranscript(aminoAcidsSequence, transcriptData));
                }
            });
            return result;
        }
    }

    private class GetSpecificTranscript implements TranscriptRetriever
    {
        final String transcriptId;

        private GetSpecificTranscript(final String transcriptId)
        {
            this.transcriptId = transcriptId;
        }

        @Override
        public Set<ProteinTranscript> getApplicableTranscripts(final GeneData geneData, final TranscriptFilter refFilter)
        {
            TranscriptData transcriptData = mEnsemblCache.getTranscriptData(geneData.GeneId, transcriptId);
            if(transcriptData == null)
            {
                String msg = format("No transcript found. Gene: %s, transcript id: %s", geneData.GeneId, transcriptId);
                throw new IllegalArgumentException(msg);
            }
            TranscriptAminoAcids aminoAcidsSequence = mTranscriptAminoAcidsMap.get(transcriptData.TransName);
            if(!refFilter.applies(aminoAcidsSequence))
            {
                String m = format("Transcript does not match. Gene: %s, transcript: %s, ref: %s", geneData.GeneId, transcriptId, refFilter);
                throw new IllegalArgumentException(m);
            }
            return Set.of(new ProteinTranscript(aminoAcidsSequence, transcriptData));
        }
    }

    private static class VariantFactory
    {
        @NotNull
        private final GeneData mGeneData;
        @NotNull
        private final AminoAcidRange mRefRange;
        @NotNull
        private final Set<ProteinTranscript> mTranscripts;

        private VariantFactory(@NotNull final GeneData geneData, @NotNull final AminoAcidRange refRange,
                @NotNull final Set<ProteinTranscript> transcript)
        {
            this.mGeneData = geneData;
            this.mRefRange = refRange;
            this.mTranscripts = transcript;
        }

        Set<ProteinVariant> buildSingleAminoAcidVariants(AminoAcidSpecification altSpec)
        {
            return mTranscripts.stream()
                    .map(t -> new SingleAminoAcidVariant(mGeneData, t.mTranscriptData, t.mAminoAcids, altSpec))
                    .collect(Collectors.toSet());
        }

        Set<ProteinVariant> buildDuplication()
        {
            return mTranscripts.stream()
                    .map(t -> new Duplication(mGeneData, t.mTranscriptData, t.mAminoAcids, mRefRange))
                    .collect(Collectors.toSet());
        }

        Set<ProteinVariant> buildFrameshift()
        {
            return mTranscripts.stream()
                    .map(t -> new Frameshift(mGeneData, t.mTranscriptData, t.mAminoAcids, mRefRange))
                    .collect(Collectors.toSet());
        }

        Set<ProteinVariant> buildStopGained()
        {
            return mTranscripts.stream()
                    .map(t -> new StopGained(mGeneData, t.mTranscriptData, t.mAminoAcids, mRefRange))
                    .collect(Collectors.toSet());
        }

        Set<ProteinVariant> buildStartLost()
        {
            return mTranscripts.stream()
                    .map(t -> new StartLost(mGeneData, t.mTranscriptData, t.mAminoAcids, mRefRange))
                    .collect(Collectors.toSet());
        }

        Set<ProteinVariant> buildDeletion()
        {
            return mTranscripts.stream()
                    .map(t -> new Deletion(mGeneData, t.mTranscriptData, t.mAminoAcids, mRefRange))
                    .collect(Collectors.toSet());
        }

        Set<ProteinVariant> buildDeletionInsertion(String altSequence)
        {
            AminoAcidSequence altAminoAcids = AminoAcidSequence.parse(altSequence);
            return mTranscripts.stream()
                    .map(t -> new DeletionInsertion(mGeneData, t.mTranscriptData, t.mAminoAcids, mRefRange, altAminoAcids))
                    .collect(Collectors.toSet());
        }

        Set<ProteinVariant> buildInsertion(String inserted)
        {
            AminoAcidSequence altAminoAcids = AminoAcidSequence.parse(inserted);
            return mTranscripts.stream()
                    .map(t -> new Insertion(mGeneData, t.mTranscriptData, t.mAminoAcids, mRefRange, altAminoAcids))
                    .collect(Collectors.toSet());
        }
    }

    VariantParser(final EnsemblDataCache ensemblDataCache, final Map<String, TranscriptAminoAcids> transcriptAminoAcidsMap)
    {
        mEnsemblCache = ensemblDataCache;
        mTranscriptAminoAcidsMap = transcriptAminoAcidsMap;
    }

    Set<ProteinVariant> parseGeneVariants(@NotNull String gene, @NotNull String variant)
    {
        return parseGeneVariant(gene, variant, new FilterTranscripts(true));
    }

    ProteinVariant parseGeneVariant(@NotNull String gene, @NotNull String variant)
    {
        Set<ProteinVariant> variants = parseGeneVariant(gene, variant, new FilterTranscripts(false));
        Preconditions.checkArgument(variants.size() == 1);
        return variants.iterator().next();
    }

    ProteinVariant parseGeneVariant(@NotNull String gene, @NotNull String specificTranscript, @NotNull String variant)
    {
        return parseGeneVariant(gene, variant, new GetSpecificTranscript(specificTranscript)).iterator().next();
    }

    ProteinVariant parse(@NotNull String expression)
    {
        String[] geneVar = extractGeneAndVariant(expression);
        String variantDescription = geneVar[0];
        if(geneVar[1].startsWith("p."))
        {
            variantDescription = geneVar[1].substring(2);
        }
        return parseGeneVariant(geneVar[0], variantDescription);
    }

    private Set<ProteinVariant> parseGeneVariant(@NotNull String gene, @NotNull String variant, TranscriptRetriever transcriptRetriever)
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

    private VariantFactory parseFactory(@NotNull String gene, @NotNull String variantDescription,
            final TranscriptRetriever transcriptRetriever)
    {
        GeneData geneData = lookupGene(gene);
        String rangeDescription = variantDescription.substring(0, variantDescription.length() - 3);
        AminoAcidRange refRange = getAminoAcidRange(variantDescription, rangeDescription);
        Set<ProteinTranscript> transcriptData = transcriptRetriever.getApplicableTranscripts(geneData, refRange);
        return new VariantFactory(geneData, refRange, transcriptData);
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
