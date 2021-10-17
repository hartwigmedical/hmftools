package com.hartwig.hmftools.common.codon;

import static com.hartwig.hmftools.common.codon.Codons.isStopCodon;
import static com.hartwig.hmftools.common.codon.Nucleotides.DNA_BASE_LIST;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNull;

import java.util.List;
import java.util.Map;
import java.util.Set;

import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.genome.region.Strand;

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

        for(String base1 : DNA_BASE_LIST)
        {
            for(String base2 : DNA_BASE_LIST)
            {
                for(String base3 : DNA_BASE_LIST)
                {
                    String trinucleotide = base1 + base2 + base3;
                    assertEquals(true, aminoAcidExists(aminoAcidToTrinucleotideMap, trinucleotide));
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