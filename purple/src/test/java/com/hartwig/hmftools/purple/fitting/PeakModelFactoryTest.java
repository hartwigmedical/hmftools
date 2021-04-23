package com.hartwig.hmftools.purple.fitting;

import static org.junit.Assert.assertEquals;

import java.io.BufferedReader;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.Date;
import java.util.List;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.utils.Doubles;

import org.jetbrains.annotations.NotNull;
import org.junit.Ignore;
import org.junit.Test;

public class PeakModelFactoryTest
{

    @Test
    public void testOffset()
    {
        final PeakModelFactory victim = new PeakModelFactory(10, 0.05);
        assertEquals(0.02, victim.offset(0.07), 0.01);
        assertEquals(-0.02, victim.offset(0.08), 0.01);
    }

    @Test
    public void testMaxBucket()
    {
        final PeakModelFactory victim = new PeakModelFactory(10, 0.05);
        victim.modelPeakHistogram(8.18, Lists.newArrayList(WeightedPloidy.create(8.18, 18, 55)));
    }

    @Ignore
    @Test
    public void testPeakModelling()
    {
        long startTime = new Date().getTime();
        WeightedPloidyHistogram victim = new WeightedPloidyHistogram(10, 0.01);
        List<ModifiableWeightedPloidy> ploidies = readResource("ploidies.tsv");

        victim.peakPloidy(10, ploidies);

        System.out.println();
        PeakModelFactory factory = new PeakModelFactory(10, 0.05);
        factory.model(ploidies);
        System.out.println(new Date().getTime() - startTime);
    }

    @Test
    public void testPeakLikelihood()
    {
        final PeakModelFactory victim = new PeakModelFactory(10, 0.05);
        final WeightedPloidy ploidy = WeightedPloidy.create(2, 35, 50);
        assertEquals(0.06, victim.ploidyLikelihood(1.8, ploidy), 0.001);
    }

    @NotNull
    public static List<ModifiableWeightedPloidy> readResource(@NotNull final String file)
    {
        final InputStream inputStream = PeakModelFactoryTest.class.getResourceAsStream("/clonality/" + file);
        final List<String> lines = new BufferedReader(new InputStreamReader(inputStream)).lines().collect(Collectors.toList());

        final List<ModifiableWeightedPloidy> result = Lists.newArrayList();
        for(String line : lines)
        {
            String[] values = line.split("\t");

            ModifiableWeightedPloidy ploidy = ModifiableWeightedPloidy.create()
                    .setPloidy(Double.parseDouble(values[0]))
                    .setAlleleReadCount(Integer.parseInt(values[1]))
                    .setTotalReadCount(Integer.parseInt(values[2]))
                    .setWeight(1);

            if(Doubles.lessThan(Double.parseDouble(values[0]), 10))
            {
                result.add(ploidy);
            }
        }

        return result;
    }

}
