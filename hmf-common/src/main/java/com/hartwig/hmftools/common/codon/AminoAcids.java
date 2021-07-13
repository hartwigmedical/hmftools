package com.hartwig.hmftools.common.codon;

import java.util.List;
import java.util.Map;
import java.util.Set;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.genome.region.Strand;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public final class AminoAcids
{

    private static final Logger LOGGER = LogManager.getLogger(AminoAcids.class);

    private static final Map<String, Set<String>> AMINO_ACID_TO_TRINUCLEOTIDES_MAP = Maps.newHashMap();
    private static final Map<String, String> TRI_LETTER_AMINO_ACID_TO_SINGLE_LETTER = Maps.newHashMap();

    static
    {
        TRI_LETTER_AMINO_ACID_TO_SINGLE_LETTER.put("Ala", "A"); // Alanine
        TRI_LETTER_AMINO_ACID_TO_SINGLE_LETTER.put("Cys", "C"); // Cysteine
        TRI_LETTER_AMINO_ACID_TO_SINGLE_LETTER.put("Asp", "D"); // Aspartic Acid
        TRI_LETTER_AMINO_ACID_TO_SINGLE_LETTER.put("Glu", "E"); // Glutamic Acid
        TRI_LETTER_AMINO_ACID_TO_SINGLE_LETTER.put("Phe", "F"); // Phenylalanine

        AMINO_ACID_TO_TRINUCLEOTIDES_MAP.put("A", Sets.newHashSet("GCA", "GCC", "GCG", "GCT")); // Ala
        AMINO_ACID_TO_TRINUCLEOTIDES_MAP.put("C", Sets.newHashSet("TGC", "TGT")); // Cys
        AMINO_ACID_TO_TRINUCLEOTIDES_MAP.put("D", Sets.newHashSet("GAC", "GAT")); // Asp
        AMINO_ACID_TO_TRINUCLEOTIDES_MAP.put("E", Sets.newHashSet("GAA", "GAG")); // Glu
        AMINO_ACID_TO_TRINUCLEOTIDES_MAP.put("F", Sets.newHashSet("TTC", "TTT")); // Phe

        TRI_LETTER_AMINO_ACID_TO_SINGLE_LETTER.put("Gly", "G"); // Glycine
        TRI_LETTER_AMINO_ACID_TO_SINGLE_LETTER.put("His", "H"); // Histidine
        TRI_LETTER_AMINO_ACID_TO_SINGLE_LETTER.put("Ile", "I"); // Isoleucine
        TRI_LETTER_AMINO_ACID_TO_SINGLE_LETTER.put("Lys", "K"); // Lysine
        TRI_LETTER_AMINO_ACID_TO_SINGLE_LETTER.put("Leu", "L"); // Leucine

        AMINO_ACID_TO_TRINUCLEOTIDES_MAP.put("G", Sets.newHashSet("GGA", "GGC", "GGG", "GGT")); // Gly
        AMINO_ACID_TO_TRINUCLEOTIDES_MAP.put("H", Sets.newHashSet("CAC", "CAT")); // His
        AMINO_ACID_TO_TRINUCLEOTIDES_MAP.put("I", Sets.newHashSet("ATA", "ATC", "ATT")); // Ile
        AMINO_ACID_TO_TRINUCLEOTIDES_MAP.put("K", Sets.newHashSet("AAA", "AAG")); // Lys
        AMINO_ACID_TO_TRINUCLEOTIDES_MAP.put("L", Sets.newHashSet("TTA", "TTG", "CTA", "CTC", "CTG", "CTT")); // Leu

        TRI_LETTER_AMINO_ACID_TO_SINGLE_LETTER.put("Met", "M"); // Methionine
        TRI_LETTER_AMINO_ACID_TO_SINGLE_LETTER.put("Asn", "N"); // Asparagine
        TRI_LETTER_AMINO_ACID_TO_SINGLE_LETTER.put("Pro", "P"); // Proline
        TRI_LETTER_AMINO_ACID_TO_SINGLE_LETTER.put("Gln", "Q"); // Glutamine
        TRI_LETTER_AMINO_ACID_TO_SINGLE_LETTER.put("Arg", "R"); // Arginine

        AMINO_ACID_TO_TRINUCLEOTIDES_MAP.put("M", Sets.newHashSet("ATG")); // Met
        AMINO_ACID_TO_TRINUCLEOTIDES_MAP.put("N", Sets.newHashSet("AAC", "AAT")); // Asn
        AMINO_ACID_TO_TRINUCLEOTIDES_MAP.put("P", Sets.newHashSet("CCA", "CCC", "CCG", "CCT")); // Pro
        AMINO_ACID_TO_TRINUCLEOTIDES_MAP.put("Q", Sets.newHashSet("CAA", "CAG")); // Gln
        AMINO_ACID_TO_TRINUCLEOTIDES_MAP.put("R", Sets.newHashSet("AGA", "AGG", "CGA", "CGC", "CGG", "CGT")); // Arg

        TRI_LETTER_AMINO_ACID_TO_SINGLE_LETTER.put("Ser", "S"); // Serine
        TRI_LETTER_AMINO_ACID_TO_SINGLE_LETTER.put("Thr", "T"); // Threonine
        TRI_LETTER_AMINO_ACID_TO_SINGLE_LETTER.put("Val", "V"); // Valine
        TRI_LETTER_AMINO_ACID_TO_SINGLE_LETTER.put("Trp", "W"); // Tryptophan
        TRI_LETTER_AMINO_ACID_TO_SINGLE_LETTER.put("Tyr", "Y"); // Tyrosine

        AMINO_ACID_TO_TRINUCLEOTIDES_MAP.put("S", Sets.newHashSet("AGC", "AGT", "TCA", "TCC", "TCG", "TCT")); // Ser
        AMINO_ACID_TO_TRINUCLEOTIDES_MAP.put("T", Sets.newHashSet("ACA", "ACC", "ACG", "ACT")); // Thr
        AMINO_ACID_TO_TRINUCLEOTIDES_MAP.put("V", Sets.newHashSet("GTA", "GTC", "GTG", "GTT")); // Val
        AMINO_ACID_TO_TRINUCLEOTIDES_MAP.put("W", Sets.newHashSet("TGG")); // Trp
        AMINO_ACID_TO_TRINUCLEOTIDES_MAP.put("Y", Sets.newHashSet("TAC", "TAT")); // Tyr
    }

    @Nullable
    public static String findAminoAcidForCodon(@NotNull String codon)
    {
        if(codon.length() != 3)
        {
            LOGGER.warn("Cannot look up amino acids for non-codons: {}", codon);
            return null;
        }

        for(Map.Entry<String, Set<String>> entrySet : AMINO_ACID_TO_TRINUCLEOTIDES_MAP.entrySet())
        {
            if(entrySet.getValue().contains(codon))
            {
                return entrySet.getKey();
            }
        }
        return null;
    }

    @NotNull
    public static List<String> allTrinucleotidesForSameAminoAcid(@NotNull String trinucleotideToFind, @NotNull Strand strand)
    {
        if(trinucleotideToFind.length() != 3)
        {
            LOGGER.warn("Cannot look up amino acids for non-trinucleotides: {}", trinucleotideToFind);
            return Lists.newArrayList();
        }

        String aminoAcid = findAminoAcidForCodon(strand == Strand.FORWARD ? trinucleotideToFind : reverseAndFlip(trinucleotideToFind));

        List<String> allTrinucleotides = Lists.newArrayList();
        Set<String> trinucleotides = AMINO_ACID_TO_TRINUCLEOTIDES_MAP.get(aminoAcid);
        if(trinucleotides != null)
        {
            for(String trinucleotide : trinucleotides)
            {
                allTrinucleotides.add(strand == Strand.FORWARD ? trinucleotide : reverseAndFlip(trinucleotide));
            }
        }
        else
        {
            LOGGER.warn("Could not find amino acid for trinucleotide '{}' on {} strand", trinucleotideToFind, strand);
        }

        return allTrinucleotides;
    }

    @NotNull
    public static String forceSingleLetterProteinAnnotation(@NotNull String proteinAnnotation)
    {
        String convertedProteinAnnotation = proteinAnnotation;
        for(Map.Entry<String, String> mapping : TRI_LETTER_AMINO_ACID_TO_SINGLE_LETTER.entrySet())
        {
            convertedProteinAnnotation = convertedProteinAnnotation.replaceAll(mapping.getKey(), mapping.getValue());
        }
        return convertedProteinAnnotation;
    }

    @NotNull
    @VisibleForTesting
    static Map<String, Set<String>> aminoAcidToTrinucleotidesMap()
    {
        return AMINO_ACID_TO_TRINUCLEOTIDES_MAP;
    }

    @NotNull
    public static String reverseAndFlip(@NotNull String string)
    {
        StringBuilder stringBuilder = new StringBuilder();
        for(int i = string.length() - 1; i >= 0; i--)
        {
            stringBuilder.append(flipBase(string.charAt(i)));
        }
        return stringBuilder.toString();
    }

    private static char flipBase(char base)
    {
        switch(base)
        {
            case 'A':
                return 'T';
            case 'T':
                return 'A';
            case 'G':
                return 'C';
            case 'C':
                return 'G';
        }

        LOGGER.warn("Cannot flip invalid base '{};'", base);
        return base;
    }
}
