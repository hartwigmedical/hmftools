package com.hartwig.hmftools.pave.transval;

import static java.lang.String.format;

import static com.hartwig.hmftools.pave.transval.Checks.HGVS_FORMAT_REQUIRED;

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
        final TranscriptAminoAcids aminoAcids;
        @NotNull
        final TranscriptData transcriptData;

        private ProteinTranscript(@NotNull final TranscriptAminoAcids aminoAcids, @NotNull final TranscriptData transcriptData)
        {
            this.aminoAcids = aminoAcids;
            this.transcriptData = transcriptData;
        }
    }

    private static class VariantFactory
    {
        @NotNull
        private final GeneData geneData;
        @NotNull
        private final AminoAcidRange refRange;
        @NotNull
        private final Set<ProteinTranscript> transcriptd;

        private VariantFactory(@NotNull final GeneData geneData, @NotNull final AminoAcidRange refRange,
                @NotNull final Set<ProteinTranscript> transcriptd)
        {
            this.geneData = geneData;
            this.refRange = refRange;
            this.transcriptd = transcriptd;
        }

        Set<ProteinVariant> buildSingleAminoAcidVariants(AminoAcidSpecification altSpec)
        {
            return transcriptd.stream()
                    .map(t -> new SingleAminoAcidVariant(geneData, t.transcriptData, t.aminoAcids, altSpec))
                    .collect(Collectors.toSet());
        }

        Set<ProteinVariant> buildDuplication()
        {
            return transcriptd.stream()
                    .map(t -> new Duplication(geneData, t.transcriptData, t.aminoAcids, refRange))
                    .collect(Collectors.toSet());
        }

        Set<ProteinVariant> buildFrameshift()
        {
            return transcriptd.stream()
                    .map(t -> new Frameshift(geneData, t.transcriptData, t.aminoAcids, refRange))
                    .collect(Collectors.toSet());
        }

        Set<ProteinVariant> buildStopGained()
        {
            return transcriptd.stream()
                    .map(t -> new StopGained(geneData, t.transcriptData, t.aminoAcids, refRange))
                    .collect(Collectors.toSet());
        }

        Set<ProteinVariant> buildStartLost()
        {
            return transcriptd.stream()
                    .map(t -> new StartLost(geneData, t.transcriptData, t.aminoAcids, refRange))
                    .collect(Collectors.toSet());
        }

        Set<ProteinVariant> buildDeletion()
        {
            return transcriptd.stream()
                    .map(t -> new Deletion(geneData, t.transcriptData, t.aminoAcids, refRange))
                    .collect(Collectors.toSet());
        }

        Set<ProteinVariant> buildDeletionInsertion(String altSequence)
        {
            AminoAcidSequence altAminoAcids = AminoAcidSequence.parse(altSequence);
            return transcriptd.stream()
                    .map(t -> new DeletionInsertion(geneData, t.transcriptData, t.aminoAcids, refRange, altAminoAcids))
                    .collect(Collectors.toSet());
        }

        Set<ProteinVariant> buildInsertion(String inserted)
        {
            AminoAcidSequence altAminoAcids = AminoAcidSequence.parse(inserted);
            return transcriptd.stream()
                    .map(t -> new Insertion(geneData, t.transcriptData, t.aminoAcids, refRange, altAminoAcids))
                    .collect(Collectors.toSet());
        }
    }

    public VariantParser(final EnsemblDataCache ensemblDataCache, final Map<String, TranscriptAminoAcids> transcriptAminoAcidsMap)
    {
        mEnsemblCache = ensemblDataCache;
        mTranscriptAminoAcidsMap = transcriptAminoAcidsMap;
    }

    public Set<ProteinVariant> parseGeneVariants(@NotNull String gene, @NotNull String variant)
    {
        return parseGeneVariant(gene, variant, true);
    }

    public ProteinVariant parseGeneVariant(@NotNull String gene, @NotNull String variant)
    {
        Set<ProteinVariant> variants = parseGeneVariant(gene, variant, false);
        Preconditions.checkArgument(variants.size() == 1);
        return variants.iterator().next();
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

    private Set<ProteinVariant> parseGeneVariant(@NotNull String gene, @NotNull String variant,
            boolean returnVariantForEachMatchingTranscript)
    {
        if(variant.contains("delins"))
        {
            String[] refAltParts = variant.split("delins");
            return parseFactory(gene, variant, returnVariantForEachMatchingTranscript).buildDeletionInsertion(refAltParts[1]);
        }
        if(variant.endsWith("del"))
        {
            return parseFactory(gene, variant, returnVariantForEachMatchingTranscript).buildDeletion();
        }
        if(variant.endsWith("dup"))
        {
            return parseFactory(gene, variant, returnVariantForEachMatchingTranscript).buildDuplication();
        }
        if(variant.contains("ins"))
        {
            String[] refAltParts = variant.split("ins");
            return parseFactory(gene, variant, returnVariantForEachMatchingTranscript).buildInsertion(refAltParts[1]);
        }
        if(variant.endsWith("fs"))
        {
            return parseFactory(gene, variant, returnVariantForEachMatchingTranscript).buildFrameshift();
        }
        if(variant.endsWith("*"))
        {
            return parseFactory(gene, variant, returnVariantForEachMatchingTranscript).buildStopGained();
        }
        if(variant.endsWith("?"))
        {
            return parseFactory(gene, variant, returnVariantForEachMatchingTranscript).buildStartLost();
        }
        Pattern variationPattern = Pattern.compile("^" + BLANK_OR_AA_GROUP + NAT + SINGLE_AA_GROUP + "$");
        final Matcher matcher = Checks.matchPattern(variationPattern, variant);
        int position = Integer.parseInt(matcher.group(3));
        AminoAcidSpecification altSpec = new AminoAcidSpecification(position, matcher.group(4));
        return parseFactory(gene, variant, returnVariantForEachMatchingTranscript).buildSingleAminoAcidVariants(altSpec);
    }

    private VariantFactory parseFactory(@NotNull String gene, @NotNull String variantDescription,
            final boolean returnVariantForEachNonCanonicalTranscript)
    {
        GeneData geneData = lookupGene(gene);
        String rangeDescription = variantDescription.substring(0, variantDescription.length() - 3);
        AminoAcidRange refRange = getAminoAcidRange(variantDescription, rangeDescription);
        Set<ProteinTranscript> transcriptData = getApplicableTranscripts(geneData, refRange, returnVariantForEachNonCanonicalTranscript);
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

    Set<ProteinTranscript> getApplicableTranscripts(final GeneData geneData,
            final TranscriptFilter refFilter,
            boolean returnVariantForEachNonCanonicalTranscript)
    {
        List<TranscriptData> allTranscripts = mEnsemblCache.getTranscripts(geneData.GeneId);
        Set<ProteinTranscript> matchingRef = filter(allTranscripts, refFilter);
        if(matchingRef.isEmpty())
        {
            String msg = format("No transcript found for gene: %s matching ref: %s", geneData.GeneId, refFilter);
            throw new IllegalArgumentException(msg);
        }
        ProteinTranscript canonical = matchingRef.stream()
                .filter(pt -> pt.transcriptData.IsCanonical)
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
