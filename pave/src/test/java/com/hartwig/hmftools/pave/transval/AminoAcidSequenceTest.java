package com.hartwig.hmftools.pave.transval;

import static org.junit.Assert.assertEquals;

import java.util.List;

import org.junit.Test;

public class AminoAcidSequenceTest extends TransvalTestBase
{
    @Test
    public void sequenceTest()
    {
        List<AminoAcid> acids = List.of(aa("A"), aa("W"), aa("R"));
        assertEquals("AWR", new AminoAcidSequence(acids).sequence());
    }

    @Test
    public void parseTest()
    {
        assertEquals("A", AminoAcidSequence.parse("A").sequence());
        assertEquals("A", AminoAcidSequence.parse("Ala").sequence());
        assertEquals("R", AminoAcidSequence.parse("Arg").sequence());
        assertEquals("RA", AminoAcidSequence.parse("ArgA").sequence());
        assertEquals("YWV", AminoAcidSequence.parse("YWV").sequence());
        assertEquals("YWV", AminoAcidSequence.parse("TyrWV").sequence());
        assertEquals("YWV", AminoAcidSequence.parse("TyrWVal").sequence());
        assertEquals("YWV", AminoAcidSequence.parse("TyrTrpVal").sequence());
    }
}
