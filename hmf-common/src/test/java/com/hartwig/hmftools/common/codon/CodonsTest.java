package com.hartwig.hmftools.common.codon;

import static org.junit.Assert.assertEquals;

import org.junit.Test;

public class CodonsTest {

    @Test
    public void testAsCodon() {
        assertEquals('Y', Codons.aminoAcid("TAT"));
    }

    @Test
    public void testStopCodon() {
        assertEquals('X', Codons.aminoAcid("TAA"));
        assertEquals('X', Codons.aminoAcid("TAG"));
        assertEquals('X', Codons.aminoAcid("TGA"));
    }

    @Test
    public void testAsCodonString() {
        assertEquals("Y", Codons.aminoAcids("TAT"));
        assertEquals("Y", Codons.aminoAcids("TATG"));
        assertEquals("Y", Codons.aminoAcids("TATGA"));
        assertEquals("YD", Codons.aminoAcids("TATGAT"));
    }

    @Test
    public void testAminoAcidToCodon() {
        assertEquals("TAT", Codons.codon('Y'));
        assertEquals("TAT", Codons.codon("Y"));
        assertEquals("TATTGA", Codons.codon("YX"));
    }

}
