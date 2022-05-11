package com.hartwig.hmftools.common.codon;

import static com.hartwig.hmftools.common.codon.Codons.START_CODON;
import static com.hartwig.hmftools.common.codon.Codons.STOP_CODON_1;
import static com.hartwig.hmftools.common.codon.Codons.STOP_CODON_2;
import static com.hartwig.hmftools.common.codon.Codons.STOP_CODON_3;
import static com.hartwig.hmftools.common.codon.Codons.UNKNOWN;
import static com.hartwig.hmftools.common.codon.Codons.codonToAminoAcid;
import static com.hartwig.hmftools.common.codon.Codons.isStopCodon;

import java.util.List;
import java.util.Map;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Lists;
import com.google.common.collect.Maps;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public final class AminoAcids
{
    // allow amino acids back to possible codons
    public static final Map<String,List<String>> AMINO_ACID_TO_CODON_MAP = Maps.newHashMap();

    // long amino-acid name to single letter
    public static final Map<String, String> TRI_LETTER_AMINO_ACID_TO_SINGLE_LETTER = Maps.newHashMap();

    static
    {
        TRI_LETTER_AMINO_ACID_TO_SINGLE_LETTER.put("Ala", "A"); // Alanine
        TRI_LETTER_AMINO_ACID_TO_SINGLE_LETTER.put("Cys", "C"); // Cysteine
        TRI_LETTER_AMINO_ACID_TO_SINGLE_LETTER.put("Asp", "D"); // Aspartic Acid
        TRI_LETTER_AMINO_ACID_TO_SINGLE_LETTER.put("Glu", "E"); // Glutamic Acid
        TRI_LETTER_AMINO_ACID_TO_SINGLE_LETTER.put("Phe", "F"); // Phenylalanine

        AMINO_ACID_TO_CODON_MAP.put("A", Lists.newArrayList("GCT", "GCA", "GCC", "GCG")); // Ala
        AMINO_ACID_TO_CODON_MAP.put("C", Lists.newArrayList("TGT", "TGC")); // Cys
        AMINO_ACID_TO_CODON_MAP.put("D", Lists.newArrayList("GAC", "GAT")); // Asp
        AMINO_ACID_TO_CODON_MAP.put("E", Lists.newArrayList("GAG", "GAA")); // Glu
        AMINO_ACID_TO_CODON_MAP.put("F", Lists.newArrayList("TTC", "TTT")); // Phe

        TRI_LETTER_AMINO_ACID_TO_SINGLE_LETTER.put("Gly", "G"); // Glycine
        TRI_LETTER_AMINO_ACID_TO_SINGLE_LETTER.put("His", "H"); // Histidine
        TRI_LETTER_AMINO_ACID_TO_SINGLE_LETTER.put("Ile", "I"); // Isoleucine
        TRI_LETTER_AMINO_ACID_TO_SINGLE_LETTER.put("Lys", "K"); // Lysine
        TRI_LETTER_AMINO_ACID_TO_SINGLE_LETTER.put("Leu", "L"); // Leucine

        AMINO_ACID_TO_CODON_MAP.put("G", Lists.newArrayList("GGA", "GGC", "GGG", "GGT")); // Gly
        AMINO_ACID_TO_CODON_MAP.put("H", Lists.newArrayList("CAC", "CAT")); // His
        AMINO_ACID_TO_CODON_MAP.put("I", Lists.newArrayList("ATA", "ATC", "ATT")); // Ile
        AMINO_ACID_TO_CODON_MAP.put("K", Lists.newArrayList("AAA", "AAG")); // Lys
        AMINO_ACID_TO_CODON_MAP.put("L", Lists.newArrayList("TTA", "TTG", "CTA", "CTC", "CTG", "CTT")); // Leu

        TRI_LETTER_AMINO_ACID_TO_SINGLE_LETTER.put("Met", "M"); // Methionine
        TRI_LETTER_AMINO_ACID_TO_SINGLE_LETTER.put("Asn", "N"); // Asparagine
        TRI_LETTER_AMINO_ACID_TO_SINGLE_LETTER.put("Pro", "P"); // Proline
        TRI_LETTER_AMINO_ACID_TO_SINGLE_LETTER.put("Gln", "Q"); // Glutamine
        TRI_LETTER_AMINO_ACID_TO_SINGLE_LETTER.put("Arg", "R"); // Arginine

        AMINO_ACID_TO_CODON_MAP.put("M", Lists.newArrayList(START_CODON)); // Met
        AMINO_ACID_TO_CODON_MAP.put("N", Lists.newArrayList("AAC", "AAT")); // Asn
        AMINO_ACID_TO_CODON_MAP.put("P", Lists.newArrayList("CCA", "CCC", "CCG", "CCT")); // Pro
        AMINO_ACID_TO_CODON_MAP.put("Q", Lists.newArrayList("CAA", "CAG")); // Gln
        AMINO_ACID_TO_CODON_MAP.put("R", Lists.newArrayList("AGG", "AGA", "CGA", "CGC", "CGG", "CGT")); // Arg

        TRI_LETTER_AMINO_ACID_TO_SINGLE_LETTER.put("Ser", "S"); // Serine
        TRI_LETTER_AMINO_ACID_TO_SINGLE_LETTER.put("Thr", "T"); // Threonine
        TRI_LETTER_AMINO_ACID_TO_SINGLE_LETTER.put("Val", "V"); // Valine
        TRI_LETTER_AMINO_ACID_TO_SINGLE_LETTER.put("Trp", "W"); // Tryptophan
        TRI_LETTER_AMINO_ACID_TO_SINGLE_LETTER.put("Tyr", "Y"); // Tyrosine

        AMINO_ACID_TO_CODON_MAP.put("S", Lists.newArrayList("TCA", "AGC", "AGT", "TCC", "TCG", "TCT")); // Ser
        AMINO_ACID_TO_CODON_MAP.put("T", Lists.newArrayList("ACA", "ACC", "ACG", "ACT")); // Thr
        AMINO_ACID_TO_CODON_MAP.put("V", Lists.newArrayList("GTG", "GTA", "GTC", "GTT")); // Val
        AMINO_ACID_TO_CODON_MAP.put("W", Lists.newArrayList("TGG")); // Trp
        AMINO_ACID_TO_CODON_MAP.put("Y", Lists.newArrayList("TAC", "TAT")); // Tyr

        AMINO_ACID_TO_CODON_MAP.put("X", Lists.newArrayList(STOP_CODON_1, STOP_CODON_2, STOP_CODON_3)); // Tyr
    }

    private static final Logger LOGGER = LogManager.getLogger(AminoAcids.class);

    @Nullable
    public static String findAminoAcidForCodon(@NotNull String codon)
    {
        // only diff is ignores the stop codon
        if(codon.length() != 3)
        {
            LOGGER.warn("Cannot look up amino acids for non-codons: {}", codon);
            return null;
        }

        char aminoAcid = codonToAminoAcid(codon);

        if(aminoAcid == UNKNOWN || isStopCodon(codon))
            return null;

        return String.valueOf(aminoAcid);
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
    static Map<String,List<String>> aminoAcidToTrinucleotidesMap()
    {
        return AMINO_ACID_TO_CODON_MAP;
    }
}
