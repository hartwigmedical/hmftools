package com.hartwig.hmftools.common.codon;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNull;

import java.util.Map;
import java.util.Set;

import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.genome.region.Strand;

import org.jetbrains.annotations.NotNull;
import org.junit.Test;

public class AminoAcidFunctionsTest {

    private static final Set<String> STOP_CODONS = Sets.newHashSet("TGA", "TAA", "TAG");

    @Test
    public void canLookupAminoAcidForTrinucleotide() {
        assertEquals("E", AminoAcids.findAminoAcidForCodon("GAA"));
        assertNull(AminoAcids.findAminoAcidForCodon("TGA"));
    }

    @Test
    public void allAminoAcidsAreFoundAndAreUnique() {
        Map<String, Set<String>> aminoAcidToTrinucleotideMap = AminoAcids.aminoAcidToTrinucleotidesMap();
        Set<String> bases = Sets.newHashSet("A", "T", "G", "C");

        for (String base1 : bases) {
            for (String base2 : bases) {
                for (String base3 : bases) {
                    String trinucleotide = base1 + base2 + base3;
                    boolean expectsToExist = !STOP_CODONS.contains(trinucleotide);
                    assertEquals(expectsToExist, aminoAcidExists(aminoAcidToTrinucleotideMap, trinucleotide));
                }
            }
        }
    }

    private static boolean aminoAcidExists(@NotNull Map<String, Set<String>> aminoAcidToTrinucleotideMap, @NotNull String trinucleotide) {
        for (Map.Entry<String, Set<String>> entry : aminoAcidToTrinucleotideMap.entrySet()) {
            if (entry.getValue().contains(trinucleotide)) {
                return true;
            }
        }
        return false;
    }

    @Test
    public void canLookupTrinucleotides() {
        // Serine (S)
        assertEquals(6, AminoAcids.allTrinucleotidesForSameAminoAcid("TCT", Strand.FORWARD).size());
        assertEquals(6, AminoAcids.allTrinucleotidesForSameAminoAcid("AGA", Strand.REVERSE).size());

        // Valine (V)
        assertEquals(4, AminoAcids.allTrinucleotidesForSameAminoAcid("GTC", Strand.FORWARD).size());
        assertEquals(4, AminoAcids.allTrinucleotidesForSameAminoAcid("GAC", Strand.REVERSE).size());

        // Tyrosine (Y)
        assertEquals(2, AminoAcids.allTrinucleotidesForSameAminoAcid("TAC", Strand.FORWARD).size());
        assertEquals(2, AminoAcids.allTrinucleotidesForSameAminoAcid("GTA", Strand.REVERSE).size());

        // Does not exist -> no trinucleotides found!
        assertEquals(0, AminoAcids.allTrinucleotidesForSameAminoAcid("???", Strand.FORWARD).size());
        assertEquals(0, AminoAcids.allTrinucleotidesForSameAminoAcid("???", Strand.REVERSE).size());

        // No trinucleotide -> return none.
        assertEquals(0, AminoAcids.allTrinucleotidesForSameAminoAcid("TCTC", Strand.FORWARD).size());
        assertEquals(0, AminoAcids.allTrinucleotidesForSameAminoAcid("GAGA", Strand.REVERSE).size());
    }

    @Test
    public void canForceSingleLetterProteinAnnotation() {
        assertEquals("p.N334K", AminoAcids.forceSingleLetterProteinAnnotation("p.Asn334Lys"));
    }
}