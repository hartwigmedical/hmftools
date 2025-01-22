package com.hartwig.hmftools.pave.transval;

import static org.junit.Assert.assertEquals;

import org.junit.Test;

public class DeletionInsertionTest extends TransvalTestBase
{

    @Test
    public void referenceAminoAcidsTest()
    {
        assertEquals("M", this.di("ADCK2:p.M1_M1delinsKQ").referenceAminoAcids());
        assertEquals("MVAP", this.di("ADCK2:p.M1_P4delinsK").referenceAminoAcids());
        assertEquals("EAT", this.di("ADCK2:p.E301_T303delinsQQ").referenceAminoAcids());
        assertEquals("PX", this.di("ADCK2:p.E626_T627delinsQQ").referenceAminoAcids());
    }

    @Test
    public void altAminoAcidsTest()
    {
        assertEquals("KQ", this.di("ADCK2:p.M1_M1delinsKQ").altAminoAcidSequence());
        assertEquals("K", this.di("ADCK2:p.M1_P4delinsK").altAminoAcidSequence());
    }

    @Test
    public void positionTest()
    {
        assertEquals(1, this.di("ADCK2:p.M1_M1delinsKQ").positionOfFirstAlteredCodon());
        assertEquals(1, this.di("ADCK2:p.M1_M1delinsKQ").positionOfLastAlteredCodon());
        assertEquals(111, this.di("VHL:p.S111_L116delinsKQ").positionOfFirstAlteredCodon());
        assertEquals(116, this.di("VHL:p.S111_L116delinsKQ").positionOfLastAlteredCodon());
    }

    @Test
    public void referenceBasesTest()
    {
        assertEquals(seq("ATG", null), this.di("VHL:p.M1_M1delinsKQ").referenceBases(genome));
        assertEquals(seq("ATGCCC", null), this.di("VHL:p.M1_P2delinsKQ").referenceBases(genome));
        assertEquals(seq("ATGCCCCGG", null), this.di("VHL:p.M1_R3delinsKQ").referenceBases(genome));
        assertEquals(seq("CCCCGG", null), this.di("VHL:p.P2_R3delinsKQ").referenceBases(genome));
        assertEquals(seq("GAGAACTGGGAC", null), this.di("VHL:p.E6_D9delinsKQ").referenceBases(genome));

        assertEquals(seq("AGCTACCGA", null), this.di("VHL:p.S111_R113delinsKQ").referenceBases(genome));
        assertEquals(seq("AGCTACCGAG", "GT"), this.di("VHL:p.S111_G114delinsKQ").referenceBases(genome));
        assertEquals(seq("AGCTACCGAG", "GTCAC"), this.di("VHL:p.S111_H115delinsKQ").referenceBases(genome));
        assertEquals(seq("AGCTACCGAG", "GTCACCTT"), this.di("VHL:p.S111_L116delinsKQ").referenceBases(genome));
        assertEquals(seq("TACCGAG", "GTCACCTT"), this.di("VHL:p.Y112_L116delinsKQ").referenceBases(genome));
        assertEquals(seq("CGAG", "GTCACCTT"), this.di("VHL:p.R113_L116delinsKQ").referenceBases(genome));
        assertEquals(seq("G", "GTCACCTT"), this.di("VHL:p.G114_L116delinsKQ").referenceBases(genome));

        assertEquals(seq("CACTGTGTCCCCGACTAC", null), this.di("ZYX:p.H491_Y496delinsW").referenceBases(genome));
        assertEquals(seq("TGTGTCCCCGACTACCAC", null), this.di("ZYX:p.C492_H497delinsW").referenceBases(genome));
        assertEquals(seq("GTCCCCGACTACCACAA", "G"), this.di("ZYX:p.V493_K498delinsW").referenceBases(genome));
        assertEquals(seq("CCCGACTACCACAA", "GCAG"), this.di("ZYX:p.P494_Q499delinsW").referenceBases(genome));
        assertEquals(seq("GACTACCACAA", "GCAGTAC"), this.di("ZYX:p.D495_Y500delinsW").referenceBases(genome));
    }

    private DeletionInsertion di(String definition)
    {
        return transval.variationParser().parseDeletionInsertion(definition);
    }

    private SplitSequence seq(String left, String right)
    {
        return new SplitSequence(left, right);
    }
}
