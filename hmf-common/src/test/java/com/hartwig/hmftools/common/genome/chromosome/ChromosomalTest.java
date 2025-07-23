package com.hartwig.hmftools.common.genome.chromosome;

import static junit.framework.TestCase.assertEquals;

import org.jetbrains.annotations.NotNull;
import org.junit.Test;

public class ChromosomalTest
{
    @Test
    public void chrTest()
    {
        assertEquals(HumanChromosome._1, new IntChr(1).chr());
    }

    @Test(expected = IllegalArgumentException.class)
    public void testUnknownChromosome()
    {

        new IntChr(123).chr();
    }
}

record IntChr(int value) implements Chromosomal
{
    @NotNull
    @Override
    public String chromosome()
    {
        return "" + value;
    }
}
