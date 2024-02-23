package com.hartwig.hmftools.lilac.fragment;

import static junit.framework.TestCase.assertEquals;

import com.google.common.collect.Lists;

import org.junit.Test;

public class NucleotideFragmentFactoryTest
{
    @Test
    public void testCreateNucleotidesFromAminoAcid()
    {
        assertEquals(Lists.newArrayList("T", "A", "A"), NucleotideFragmentFactory.createNucleotidesFromAminoAcid("X"));
        assertEquals(Lists.newArrayList(".", ".", "."), NucleotideFragmentFactory.createNucleotidesFromAminoAcid("."));
        assertEquals(Lists.newArrayList("A", "G", "CTAA"), NucleotideFragmentFactory.createNucleotidesFromAminoAcid("SX"));
    }
}
