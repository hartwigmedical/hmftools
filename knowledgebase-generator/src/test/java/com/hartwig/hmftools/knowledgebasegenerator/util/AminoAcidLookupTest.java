package com.hartwig.hmftools.knowledgebasegenerator.util;

import static org.junit.Assert.assertEquals;

import java.util.Map;
import java.util.Set;

import com.google.common.collect.Sets;

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

    private int aminoAcidCount(@NotNull Map<String, Set<String>> aminoAcidToTrinucleotideMap, @NotNull String trinucleotide) {
        int count = 0;
        for (Map.Entry<String, Set<String>> entry : aminoAcidToTrinucleotideMap.entrySet()) {
            if (entry.getValue().contains(trinucleotide)) {
                count++;
            }
        }
        return count;
    }
}