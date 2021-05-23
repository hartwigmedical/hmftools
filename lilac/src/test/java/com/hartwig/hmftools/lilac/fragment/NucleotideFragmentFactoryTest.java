package com.hartwig.hmftools.lilac.fragment;

import static junit.framework.TestCase.assertEquals;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.lilac.hla.HlaAllele;
import com.hartwig.hmftools.lilac.seq.HlaSequence;

import org.junit.Test;

import htsjdk.samtools.SAMRecord;

public class NucleotideFragmentFactoryTest
{
    /*
    private HlaSequence reference = new HlaSequence(HlaAllele.fromString("A*01:01"), "GATTTACA..CAT..ATC");
    private HlaSequence delete = new HlaSequence(HlaAllele.fromString("A*01:02"), "--...---..---..---");
    private HlaSequence insert = new HlaSequence(HlaAllele.fromString("A*01:03"), "--------AA---TT---");
    */

    @Test
    public void testCreateNucleotidesFromAminoAcid()
    {
        assertEquals(Lists.newArrayList("T", "G", "A"), NucleotideFragmentFactory.createNucleotidesFromAminoAcid("X"));
        assertEquals(Lists.newArrayList(".", ".", "."), NucleotideFragmentFactory.createNucleotidesFromAminoAcid("."));
        assertEquals(Lists.newArrayList("A", "G", "TTGA"), NucleotideFragmentFactory.createNucleotidesFromAminoAcid("SX"));
    }
}
