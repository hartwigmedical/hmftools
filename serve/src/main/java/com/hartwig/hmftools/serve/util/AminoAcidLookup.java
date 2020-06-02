package com.hartwig.hmftools.serve.util;

import java.util.List;
import java.util.Map;
import java.util.Set;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.genome.region.Strand;

import org.jetbrains.annotations.NotNull;

public final class AminoAcidLookup {

    private static final Map<String, Set<String>> AMINO_ACID_TO_TRINUCLEOTIDES_MAP = Maps.newHashMap();

    static {
        AMINO_ACID_TO_TRINUCLEOTIDES_MAP.put("A", Sets.newHashSet("GCA", "GCC", "GCG", "GCT")); // Ala
        AMINO_ACID_TO_TRINUCLEOTIDES_MAP.put("C", Sets.newHashSet("TGC", "TGT")); // Cys
        AMINO_ACID_TO_TRINUCLEOTIDES_MAP.put("D", Sets.newHashSet("GAC", "GAT")); // Asp
        AMINO_ACID_TO_TRINUCLEOTIDES_MAP.put("E", Sets.newHashSet("GAA", "GAG")); // Glu
        AMINO_ACID_TO_TRINUCLEOTIDES_MAP.put("F", Sets.newHashSet("TTC", "TTT")); // Phe

        AMINO_ACID_TO_TRINUCLEOTIDES_MAP.put("G", Sets.newHashSet("GGA", "GGC", "GGG", "GGT")); // Gly
        AMINO_ACID_TO_TRINUCLEOTIDES_MAP.put("H", Sets.newHashSet("CAC", "CAT")); // His
        AMINO_ACID_TO_TRINUCLEOTIDES_MAP.put("I", Sets.newHashSet("ATA", "ATC", "ATT")); // Ile
        AMINO_ACID_TO_TRINUCLEOTIDES_MAP.put("K", Sets.newHashSet("AAA", "AAG")); // Lys
        AMINO_ACID_TO_TRINUCLEOTIDES_MAP.put("L", Sets.newHashSet("TTA", "TTG", "CTA", "CTC", "CTG", "CTT")); // Leu

        AMINO_ACID_TO_TRINUCLEOTIDES_MAP.put("M", Sets.newHashSet("ATG")); // Met
        AMINO_ACID_TO_TRINUCLEOTIDES_MAP.put("N", Sets.newHashSet("AAC", "AAT")); // Asn
        AMINO_ACID_TO_TRINUCLEOTIDES_MAP.put("P", Sets.newHashSet("CCA", "CCC", "CCG", "CCT")); // Pro
        AMINO_ACID_TO_TRINUCLEOTIDES_MAP.put("Q", Sets.newHashSet("CAA", "CAG")); // Gln
        AMINO_ACID_TO_TRINUCLEOTIDES_MAP.put("R", Sets.newHashSet("AGA", "AGG", "CGA", "CGC", "CGG", "CGT")); // Arg

        AMINO_ACID_TO_TRINUCLEOTIDES_MAP.put("S", Sets.newHashSet("AGC", "AGT", "TCA", "TCC", "TCG", "TCT")); // Ser
        AMINO_ACID_TO_TRINUCLEOTIDES_MAP.put("T", Sets.newHashSet("ACA", "ACC", "ACG", "ACT")); // Thr
        AMINO_ACID_TO_TRINUCLEOTIDES_MAP.put("V", Sets.newHashSet("GTA", "GTC", "GTG", "GTT")); // Val
        AMINO_ACID_TO_TRINUCLEOTIDES_MAP.put("W", Sets.newHashSet("TGG")); // Trp
        AMINO_ACID_TO_TRINUCLEOTIDES_MAP.put("Y", Sets.newHashSet("TAC", "TAT")); // Tyr
    }

    @NotNull
    public static List<String> allTrinucleotidesForSameAminoAcid(@NotNull String trinucleotideToFind, @NotNull Strand strand) {
        assert trinucleotideToFind.length() == 3;
        String strandCorrectedTrinucleotide = strand == Strand.FORWARD ? trinucleotideToFind : reverse(trinucleotideToFind);

        List<String> allTrinucleotides = Lists.newArrayList();
        for (Map.Entry<String, Set<String>> aminoAcidEntry : AMINO_ACID_TO_TRINUCLEOTIDES_MAP.entrySet()) {
            Set<String> trinucleotides = aminoAcidEntry.getValue();
            if (trinucleotides.contains(strandCorrectedTrinucleotide)) {
                for (String trinucleotide : trinucleotides) {
                    allTrinucleotides.add(strand == Strand.FORWARD ? trinucleotide : reverse(trinucleotide));
                }
            }
        }

        return allTrinucleotides;
    }

    @NotNull
    @VisibleForTesting
    static Map<String, Set<String>> aminoAcidToTrinucleotidesMap() {
        return AMINO_ACID_TO_TRINUCLEOTIDES_MAP;
    }

    @NotNull
    private static String reverse(@NotNull String string) {
        StringBuilder stringBuilder = new StringBuilder();
        for (int i = string.length() - 1; i >= 0; i--) {
            stringBuilder.append(string.charAt(i));
        }
        return stringBuilder.toString();
    }
}
