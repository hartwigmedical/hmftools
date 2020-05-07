package com.hartwig.hmftools.serve.util;

import static org.junit.Assert.assertEquals;

import java.util.Map;
import java.util.Set;

import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.genome.region.Strand;

import org.jetbrains.annotations.NotNull;
import org.junit.Test;

public class AminoAcidLookupTest {

    private static final Set<String> STOP_CODONS = Sets.newHashSet("TGA", "TAA", "TAG");

    @Test
    public void allAminoAcidsAreFoundAndAreUnique() {
        Map<String, Set<String>> aminoAcidToTrinucleotideMap = AminoAcidLookup.aminoAcidToTrinucleotidesMap();
        Set<String> bases = Sets.newHashSet("A", "T", "G", "C");

        for (String base1 : bases) {
            for (String base2 : bases) {
                for (String base3 : bases) {
                    String trinucleotide = base1 + base2 + base3;
                    int expectedCount = STOP_CODONS.contains(trinucleotide) ? 0 : 1;
                    System.out.println(trinucleotide);
                    assertEquals(expectedCount, aminoAcidCount(aminoAcidToTrinucleotideMap, trinucleotide));
                }
            }
        }
    }

    private static int aminoAcidCount(@NotNull Map<String, Set<String>> aminoAcidToTrinucleotideMap, @NotNull String trinucleotide) {
        int count = 0;
        for (Map.Entry<String, Set<String>> entry : aminoAcidToTrinucleotideMap.entrySet()) {
            if (entry.getValue().contains(trinucleotide)) {
                count++;
            }
        }
        return count;
    }

    @Test
    public void canLookupTrinucleotides() {
        // Serine (S)
        assertEquals(6, AminoAcidLookup.allTrinucleotidesForSameAminoAcid("TCT", Strand.FORWARD).size());
        assertEquals(6, AminoAcidLookup.allTrinucleotidesForSameAminoAcid("TCT", Strand.REVERSE).size());

        // Valine (V)
        assertEquals(4, AminoAcidLookup.allTrinucleotidesForSameAminoAcid("GTC", Strand.FORWARD).size());
        assertEquals(4, AminoAcidLookup.allTrinucleotidesForSameAminoAcid("CTG", Strand.REVERSE).size());

        // Tyrosine (Y)
        assertEquals(2, AminoAcidLookup.allTrinucleotidesForSameAminoAcid("TAC", Strand.FORWARD).size());
        assertEquals(2, AminoAcidLookup.allTrinucleotidesForSameAminoAcid("CAT", Strand.REVERSE).size());

        // Does not exist -> no trinucleotides found!
        assertEquals(0, AminoAcidLookup.allTrinucleotidesForSameAminoAcid("???", Strand.FORWARD).size());
        assertEquals(0, AminoAcidLookup.allTrinucleotidesForSameAminoAcid("???", Strand.REVERSE).size());
    }
}