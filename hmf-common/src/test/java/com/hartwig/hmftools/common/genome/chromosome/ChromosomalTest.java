package com.hartwig.hmftools.common.genome.chromosome;

import static junit.framework.TestCase.assertEquals;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

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

    @Test
    public void toPerChromosomeListsTest()
    {
        assertEquals(new HashMap<>(), Chromosomal.toPerChromosomeLists(new ArrayList<>()));

        List<Chromosomal> input = new ArrayList<>();
        final ChrBase cb1_1 = new ChrBase(1, 100);
        input.add(cb1_1);
        final ChrBase cb2_1 = new ChrBase(2, 100);
        input.add(cb2_1);
        final ChrBase cb3_1 = new ChrBase(3, 100);
        input.add(cb3_1);
        final ChrBase cb1_2 = new ChrBase(1, 200);
        input.add(cb1_2);
        final ChrBase cb2_2 = new ChrBase(2, 200);
        input.add(cb2_2);
        final ChrBase cb3_2 = new ChrBase(3, 200);
        input.add(cb3_2);
        final ChrBase cb1_3 = new ChrBase(1, 300);
        input.add(cb1_3);
        final ChrBase cb2_3 = new ChrBase(2, 300);
        input.add(cb2_3);
        final ChrBase cb3_3 = new ChrBase(3, 300);
        input.add(cb3_3);
        Map<Chromosome, List<ChrBase>> expected = new HashMap<>();
        expected.put(cb1_1.chr(), new ArrayList<>(List.of(cb1_1, cb1_2, cb1_3)));
        expected.put(cb2_1.chr(), new ArrayList<>(List.of(cb2_1, cb2_2, cb2_3)));
        expected.put(cb3_1.chr(), new ArrayList<>(List.of(cb3_1, cb3_2, cb3_3)));

        assertEquals(expected, Chromosomal.toPerChromosomeLists(input));
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
record ChrBase(int value, int base) implements Chromosomal
{
    @NotNull
    @Override
    public String chromosome()
    {
        return "" + value;
    }
}
