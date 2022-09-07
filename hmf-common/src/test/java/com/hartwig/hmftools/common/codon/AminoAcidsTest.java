package com.hartwig.hmftools.common.codon;

import static com.hartwig.hmftools.common.codon.Nucleotides.DNA_BASES;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNull;
import static org.junit.Assert.assertTrue;

import java.util.List;
import java.util.Map;

import org.jetbrains.annotations.NotNull;
import org.junit.Test;

public class AminoAcidsTest
{
    @Test
    public void canLookupAminoAcidForTrinucleotide()
    {
        assertEquals("E", AminoAcids.findAminoAcidForCodon("GAA"));
        assertNull(AminoAcids.findAminoAcidForCodon("TGA"));
    }

    @Test
    public void allAminoAcidsAreFoundAndAreUnique()
    {
        Map<String, List<String>> aminoAcidToTrinucleotideMap = AminoAcids.aminoAcidToTrinucleotidesMap();

        for(char base1 : DNA_BASES)
        {
            for(char base2 : DNA_BASES)
            {
                for(char base3 : DNA_BASES)
                {
                    String trinucleotide = Character.toString(base1) + base2 + base3;
                    assertTrue(aminoAcidExists(aminoAcidToTrinucleotideMap, trinucleotide));
                }
            }
        }
    }

    private static boolean aminoAcidExists(@NotNull Map<String,List<String>> aminoAcidToTrinucleotideMap, @NotNull String trinucleotide)
    {
        for(Map.Entry<String,List<String>> entry : aminoAcidToTrinucleotideMap.entrySet())
        {
            if(entry.getValue().contains(trinucleotide))
            {
                return true;
            }
        }
        return false;
    }

    @Test
    public void canForceSingleLetterProteinAnnotation()
    {
        assertEquals("p.N334K", AminoAcids.forceSingleLetterProteinAnnotation("p.Asn334Lys"));
    }
}