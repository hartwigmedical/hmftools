package com.hartwig.hmftools.common.codon;

import static com.hartwig.hmftools.common.codon.Nucleotides.DNA_BASES;

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
        assertEquals("TAC", Codons.aminoAcidToCodon('Y'));
        assertEquals("TAC", Codons.aminoAcidsToCodons("Y"));
        assertEquals("TACTAA", Codons.aminoAcidsToCodons("YX"));
    }

    @Test
    public void testAllCodonToAminoAcid()
    {
        for(char base1 : DNA_BASES)
        {
            for(char base2 : DNA_BASES)
            {
                for(char base3 : DNA_BASES)
                {
                    String codon = "" + base1 + base2 + base3;
                    assertEquals(codonToAminoAcidCheck(codon), Codons.codonToAminoAcid(codon, 0));
                }
            }
        }
    }

    private static char codonToAminoAcidCheck(final String codon)
    {
        if(Codons.isStopCodon(codon))
            return Codons.STOP_AMINO_ACID;

        if(Codons.isStartCodon(codon))
            return Codons.START_AMINO_ACID;

        switch(codon)
        {
            // SECOND BASE T
            case "TTT":
            case "TTC":
                return 'F';
            case "TTA":
            case "TTG":
            case "CTT":
            case "CTC":
            case "CTA":
            case "CTG":
                return 'L';
            case "ATT":
            case "ATC":
            case "ATA":
                return 'I';
            case "GTT":
            case "GTC":
            case "GTA":
            case "GTG":
                return 'V';

            // SECOND BASE C
            case "TCT":
            case "TCC":
            case "TCA":
            case "TCG":
                return 'S';
            case "CCT":
            case "CCC":
            case "CCA":
            case "CCG":
                return 'P';
            case "ACT":
            case "ACC":
            case "ACA":
            case "ACG":
                return 'T';
            case "GCT":
            case "GCC":
            case "GCA":
            case "GCG":
                return 'A';

            // SECOND BASE A
            case "TAT":
            case "TAC":
                return 'Y';
            case "CAT":
            case "CAC":
                return 'H';
            case "CAA":
            case "CAG":
                return 'Q';
            case "AAT":
            case "AAC":
                return 'N';
            case "AAA":
            case "AAG":
                return 'K';
            case "GAT":
            case "GAC":
                return 'D';
            case "GAA":
            case "GAG":
                return 'E';

            // SECOND BASE G
            case "TGT":
            case "TGC":
                return 'C';
            case "TGG":
                return 'W';
            case "CGT":
            case "CGC":
            case "CGA":
            case "CGG":
                return 'R';
            case "AGT":
            case "AGC":
                return 'S';
            case "AGA":
            case "AGG":
                return 'R';
            case "GGT":
            case "GGC":
            case "GGA":
            case "GGG":
                return 'G';
        }

        return Codons.UNKNOWN;
    }
}
