package com.hartwig.hmftools.pave.transval;

import java.util.List;
import java.util.Map;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import java.util.stream.Collectors;

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
        if(expression.contains("delins"))
        {
            return parseDeletionInsertion(expression);
        }
        return parseSingleAminoAcidVariant(expression);
    }

    public ProteinVariant parseExpressionForGene(@NotNull String gene, @NotNull String expression)
    {
        if(expression.contains("delins"))
        {
            return parseDeletionInsertion(gene, expression);
        }
        return parseSingleAminoAcidVariant(gene, expression);
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

    public DeletionInsertion parseDeletionInsertion(@NotNull String gene, @NotNull String variantDescription)
    {
        GeneData geneData = lookupGene(gene);
        return buildDeletionInsertion(variantDescription, geneData);
    }

    public DeletionInsertion parseDeletionInsertion(String input)
    {
        String[] geneVar = extractGeneAndVariant(input);
        GeneData geneData = lookupGene(geneVar[0]);
        String description = extractDescription(geneVar)[1];
        return buildDeletionInsertion(description, geneData);
    }

    @NotNull
    private DeletionInsertion buildDeletionInsertion(final String description, final GeneData geneData)
    {
        TranscriptData canonicalTranscript = getCanonicalTranscriptData(geneData);
        TranscriptAminoAcids aminoAcidsSequence = lookupTranscriptAminoAcids(canonicalTranscript, false);

        Pattern variationPattern =
                Pattern.compile("^" + BLANK_OR_AA_GROUP + NAT + "_" + BLANK_OR_AA_GROUP + NAT + "delins([a-zA-Z]*)" + "$");
        final Matcher matcher = matchPattern(variationPattern, description);

        ensureAminoAcidOrBlank(matcher.group(1));
        int startPosition = Integer.parseInt(matcher.group(3));
        ensureAminoAcid(matcher.group(4));
        int endPosition = Integer.parseInt(matcher.group(6));

        String variant = matcher.group(7);
        Pattern variantListPattern = Pattern.compile("([A-Z][a-z]*)");
        final Matcher variantListMatcher = variantListPattern.matcher(variant);
        StringBuilder aaListBuilder = new StringBuilder();
        while(variantListMatcher.find())
        {
            String aa = variantListMatcher.group();
            if(isNotValidAminoAcidIdentifier(aa))
            {
                throw new IllegalArgumentException("Not a known amino acid: " + aa);
            }
            String singleLetterName = AminoAcids.forceSingleLetterProteinAnnotation(aa);
            aaListBuilder.append(singleLetterName);
        }
        String variantString = aaListBuilder.toString();

        int length = endPosition - startPosition + 1;
        if(length < 1)
        {
            throw new IllegalArgumentException("End position must not be before start position");
        }
        return new DeletionInsertion(geneData, canonicalTranscript, aminoAcidsSequence, startPosition, length, variantString);
    }

    @NotNull
    private SingleAminoAcidVariant buildSingleAminoAcidVariant(@NotNull final String variantDescription,
            final GeneData geneData)
    {
        Pattern variationPattern = Pattern.compile("^" + BLANK_OR_AA_GROUP + NAT + SINGLE_AA_GROUP + "$");
        final Matcher matcher = matchPattern(variationPattern, variantDescription);

        String ref = matcher.group(1);
        ensureAminoAcidOrBlank(ref);
        String r = (ref == null || ref.isBlank()) ? null : AminoAcids.forceSingleLetterProteinAnnotation(ref);
        int position = Integer.parseInt(matcher.group(3));
        String alt = matcher.group(4);
        ensureAminoAcid(alt);
        String a = AminoAcids.forceSingleLetterProteinAnnotation(alt);
        TranscriptData transcriptData = getApplicableTranscript(geneData, position, r, a);
        TranscriptAminoAcids aminoAcidsSequence = lookupTranscriptAminoAcids(transcriptData, false);

        return new SingleAminoAcidVariant(geneData, transcriptData, aminoAcidsSequence, position, a);
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

    @NotNull
    private static Matcher matchPattern(final Pattern variationPattern, final String description)
    {
        final Matcher matcher = variationPattern.matcher(description);
        boolean matches = matcher.find();
        if(!matches)
        {
            throw new IllegalArgumentException(HGVS_FORMAT_REQUIRED);
        }
        return matcher;
    }

    public TranscriptData getApplicableTranscript(final GeneData geneData, int position, final String ref, final String alt)
    {
        List<TranscriptData> allTranscripts = mEnsemblCache.getTranscripts(geneData.GeneId);
        List<TranscriptData> possibilities = allTranscripts
                .stream()
                .filter(transcriptData -> transcriptMatches(transcriptData, position - 1, ref, alt))
                .collect(Collectors.toList());
        if(possibilities.isEmpty())
        {
            String msg = String.format("No transcript found for gene: %s, pos: %d, ref: %s, alt: %s", geneData.GeneId, position, ref, alt);
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
                                String.format("No canonical transcript, but multiple non-canonical transcripts, found for gene: %s, pos: %d, ref: %s, alt: %s", geneData.GeneId, position, ref, alt);
                        return new IllegalStateException(msg);
                    });
        }
        return possibilities.get(0);
    }

    public boolean transcriptMatches(TranscriptData transcriptData, int position, String ref, String alt)
    {
        TranscriptAminoAcids aminoAcidsSequence = lookupTranscriptAminoAcids(transcriptData, true);
        if(aminoAcidsSequence == null)
        {
            return false;
        }
        if(aminoAcidsSequence.AminoAcids.length() <= position)
        {
//            System.out.println(transcriptData.TransName + " too short to have value at position " + position);
            return false;
        }

        if (ref != null)
        {
            String transcriptValueAtPosition = aminoAcidsSequence.AminoAcids.substring(position, position + ref.length());
            if(!ref.equals(transcriptValueAtPosition))
            {
                //            System.out.println(transcriptData.TransName + " does not match ref at position " + position);
                return false;
            }
            if(alt.equals(transcriptValueAtPosition))
            {
                //            System.out.println(transcriptData.TransName + " already matches alt at position " + position);
                return false;
            }
        }
        return true;
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
        if(isNotValidAminoAcidIdentifier(variantAminoAcid))
        {
            throw new IllegalArgumentException(variantAminoAcid + " does not represent an amino acid");
        }
    }

    private void ensureAminoAcidOrBlank(final String referenceAminoAcid)
    {
        if(referenceAminoAcid != null && !referenceAminoAcid.isBlank())
        {
            ensureAminoAcid(referenceAminoAcid);
        }
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

    private boolean isNotValidAminoAcidIdentifier(String s)
    {
        if(AminoAcids.AMINO_ACID_TO_CODON_MAP.containsKey(s))
        {
            return false;
        }
        if(AminoAcids.TRI_LETTER_AMINO_ACID_TO_SINGLE_LETTER.containsKey(s))
        {
            return false;
        }
        if("Z".equals(s) || "Glx".equals(s)) // means 'Glutamine or Glutamic Acid'
        {
            return false;
        }
        // B and Asx mean 'Asparagine or Aspartic Acid'
        return !"B".equals(s) && !"Asx".equals(s);
    }
}
