package com.hartwig.hmftools.common.codon;

import static org.junit.Assert.assertEquals;

import org.junit.Test;

public class CodonsTest
{
    @Test
    public void testAsCodon()
    {
        assertEquals('Y', Codons.codonToAminoAcid("TAT"));
    }

    @Test
    public void testStopCodon()
    {
        assertEquals('X', Codons.codonToAminoAcid("TAA"));
        assertEquals('X', Codons.codonToAminoAcid("TAG"));
        assertEquals('X', Codons.codonToAminoAcid("TGA"));
    }

    @Test
    public void testAsCodonString()
    {
        assertEquals("Y", Codons.aminoAcidFromBases("TAT"));
        assertEquals("Y", Codons.aminoAcidFromBases("TATG"));
        assertEquals("Y", Codons.aminoAcidFromBases("TATGA"));
        assertEquals("YD", Codons.aminoAcidFromBases("TATGAT"));
    }

    @Test
    public void testAminoAcidToCodon()
    {
        // there are multiple solutions here and the choice depends on the order of the DNA_BASES array
        assertEquals("TAT", Codons.aminoAcidToCodon('Y'));
        assertEquals("TAT", Codons.aminoAcidsToCodons("Y"));
        assertEquals("TATTGA", Codons.aminoAcidsToCodons("YX"));
    }
}
