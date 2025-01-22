package com.hartwig.hmftools.pave.transval;

import java.util.List;
import java.util.Map;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import com.hartwig.hmftools.common.codon.AminoAcids;
import com.hartwig.hmftools.common.ensemblcache.EnsemblDataCache;
import com.hartwig.hmftools.common.gene.GeneData;
import com.hartwig.hmftools.common.gene.TranscriptAminoAcids;
import com.hartwig.hmftools.common.gene.TranscriptData;

import org.jetbrains.annotations.NotNull;

public class VariationParser
{
    static final String HGVS_FORMAT_REQUIRED = "Required format is GENE:p.XnY where X and Y are amino acids and n is an integer.";
    private final EnsemblDataCache mEnsemblCache;
    private final Map<String, TranscriptAminoAcids> mTranscriptAminoAcidsMap;
    final String SINGLE_AA = "[A-Z]([a-z][a-z])?+";
    final String SINGLE_AA_GROUP = "(" + SINGLE_AA + ")";
    final String BLANK_OR_AA_GROUP = "(" + SINGLE_AA + ")?+";
    final String NAT = "(\\d+)";


    public VariationParser(final EnsemblDataCache ensemblDataCache, final Map<String, TranscriptAminoAcids> transcriptAminoAcidsMap)
    {
        mEnsemblCache = ensemblDataCache;
        mTranscriptAminoAcidsMap = transcriptAminoAcidsMap;
    }

    public ProteinVariant parse(@NotNull String expression)
    {
        if (expression.contains("delins"))
        {
            return parseDeletionInsertion(expression);
        }
        return parseSingleAminoAcidVariant(expression);
    }

    public SingleAminoAcidVariant parseSingleAminoAcidVariant(String input)
    {
        String[] geneVar = extractGeneAndVariant(input);
        GeneData geneData = lookupGene(geneVar[0]);
        TranscriptData canonicalTranscript = getCanonicalTranscriptData(geneData);
        TranscriptAminoAcids aminoAcidsSequence = lookupTranscriptAminoAcids(canonicalTranscript);

        String[] descriptionParts = extractDescription(geneVar);
        Pattern variationPattern = Pattern.compile("^" + BLANK_OR_AA_GROUP + NAT + SINGLE_AA_GROUP + "$");
        final Matcher matcher = matchPattern(variationPattern, descriptionParts[1]);

        String referenceAminoAcid = matcher.group(1);
        ensureAminoAcidOrBlank(referenceAminoAcid);
        int position = Integer.parseInt(matcher.group(3));
        String variantAminoAcid = matcher.group(4);
        ensureAminoAcid(variantAminoAcid);
        String singleLetterName = AminoAcids.forceSingleLetterProteinAnnotation(variantAminoAcid);

        return new SingleAminoAcidVariant(geneData, canonicalTranscript, aminoAcidsSequence, position, singleLetterName);
    }

    public DeletionInsertion parseDeletionInsertion(String input)
    {
        String[] geneVar = extractGeneAndVariant(input);
        GeneData geneData = lookupGene(geneVar[0]);
        TranscriptData canonicalTranscript = getCanonicalTranscriptData(geneData);
        TranscriptAminoAcids aminoAcidsSequence = lookupTranscriptAminoAcids(canonicalTranscript);
        String[] descriptionParts = extractDescription(geneVar);

        Pattern variationPattern = Pattern.compile("^" + BLANK_OR_AA_GROUP + NAT + "_" + BLANK_OR_AA_GROUP + NAT + "delins([a-zA-Z]*)" + "$");
        final Matcher matcher = matchPattern(variationPattern, descriptionParts[1]);

        ensureAminoAcidOrBlank(matcher.group(1));
        int startPosition = Integer.parseInt(matcher.group(3));
        ensureAminoAcid(matcher.group(4));
        int endPosition = Integer.parseInt(matcher.group(6));

        String variant = matcher.group(7);
        Pattern variantListPattern = Pattern.compile("([A-Z][a-z]*)");
        final Matcher variantListMatcher = variantListPattern.matcher(variant);
        StringBuilder aaListBuilder = new StringBuilder();
        while (variantListMatcher.find())
        {
            String aa = variantListMatcher.group();
            if (isNotValidAminoAcidIdentifier(aa))
            {
                throw new IllegalArgumentException("Not a known amino acid: " + aa);
            }
            String singleLetterName = AminoAcids.forceSingleLetterProteinAnnotation(aa);
            aaListBuilder.append(singleLetterName);
        }
        String variantString = aaListBuilder.toString();

        int length = endPosition - startPosition + 1;
        if (length < 1)
        {
            throw new IllegalArgumentException("End position must not be before start position");
        }
        return new DeletionInsertion(geneData, canonicalTranscript, aminoAcidsSequence, startPosition, length, variantString);
    }

    @NotNull
    private TranscriptAminoAcids lookupTranscriptAminoAcids(final TranscriptData canonicalTranscript)
    {
        TranscriptAminoAcids aminoAcidsSequence = mTranscriptAminoAcidsMap.get(canonicalTranscript.TransName);
        if (aminoAcidsSequence == null)
        {
            throw new IllegalStateException("No amino acid sequence found for " + canonicalTranscript.TransName);
        }
        return aminoAcidsSequence;
    }

    @NotNull
    private static Matcher matchPattern(final Pattern variationPattern, final String description)
    {
        final Matcher matcher = variationPattern.matcher(description);
        boolean matches = matcher.find();
        if (!matches)
        {
            throw new IllegalArgumentException(HGVS_FORMAT_REQUIRED);
        }
        return matcher;
    }

    private TranscriptData getCanonicalTranscriptData(final GeneData geneData)
    {
        List<TranscriptData> transcripts = mEnsemblCache.getTranscripts(geneData.GeneId);
        return transcripts.stream()
                .filter(transcriptData -> transcriptData.IsCanonical)
                .findFirst()
                .orElseThrow(() -> new IllegalStateException("No canonical transcript found for " + geneData.GeneId));
    }

    private void ensureAminoAcid(final String variantAminoAcid)
    {
        if (isNotValidAminoAcidIdentifier(variantAminoAcid))
        {
            throw new IllegalArgumentException(variantAminoAcid + " does not represent an amino acid");
        }
    }

    private void ensureAminoAcidOrBlank(final String referenceAminoAcid)
    {
        if (referenceAminoAcid != null && !referenceAminoAcid.isBlank())
        {
            ensureAminoAcid(referenceAminoAcid);
        }
    }

    @NotNull
    private GeneData lookupGene(final String reference)
    {
        GeneData geneData = mEnsemblCache.getGeneDataByName(reference);
        if (geneData == null)
        {
            geneData = mEnsemblCache.getGeneDataById(reference);
            if (geneData == null)
            {
                throw new IllegalArgumentException(reference + " is not a known gene");
            }
        }
        return geneData;
    }

    private static String[] extractDescription(final String[] geneVar)
    {
        String[] descriptionParts = geneVar[1].split("\\.");
        if (descriptionParts.length != 2)
        {
            throw new IllegalArgumentException(HGVS_FORMAT_REQUIRED);
        }
        if (!descriptionParts[0].equals("p"))
        {
            throw new IllegalArgumentException(HGVS_FORMAT_REQUIRED);
        }
        return descriptionParts;
    }

    private static String[] extractGeneAndVariant(final String input)
    {
        String[] geneVar = input.split(":");
        if (geneVar.length != 2)
        {
            throw new IllegalArgumentException(HGVS_FORMAT_REQUIRED);
        }
        return geneVar;
    }

    private boolean isNotValidAminoAcidIdentifier(String s)
    {
        if (AminoAcids.AMINO_ACID_TO_CODON_MAP.containsKey(s))
        {
            return false;
        }
        if (AminoAcids.TRI_LETTER_AMINO_ACID_TO_SINGLE_LETTER.containsKey(s))
        {
            return false;
        }
        if ("Z".equals(s) || "Glx".equals(s)) // means 'Glutamine or Glutamic Acid'
        {
            return false;
        }
        // B and Asx mean 'Asparagine or Aspartic Acid'
        return !"B".equals(s) && !"Asx".equals(s);
    }
}
