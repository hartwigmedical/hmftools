package com.hartwig.hmftools.purple.fitting;

import static com.hartwig.hmftools.common.utils.file.FileDelimiters.TSV_DELIM;

import static org.junit.Assert.assertEquals;

import java.io.BufferedReader;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.List;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.utils.Doubles;
import com.hartwig.hmftools.purple.fittingsnv.WeightedPloidy;
import com.hartwig.hmftools.purple.fittingsnv.WeightedPloidyHistogram;

import org.jetbrains.annotations.NotNull;
import org.junit.Test;

public class WeightedPloidyHistogramTest
{
    @Test
    public void testMaxBucket()
    {
        final double max = 5;
        final WeightedPloidyHistogram victim = new WeightedPloidyHistogram(max, 0.01);
        victim.histogram(Lists.newArrayList(WeightedPloidy.create(max, 20, 40)));
    }

    @Test
    public void testBucketWithOffset()
    {
        final WeightedPloidyHistogram negativeOffset = new WeightedPloidyHistogram(10, 0.05, -0.02);
        assertEquals(1, negativeOffset.bucket(0.05));
        assertEquals(2, negativeOffset.bucket(0.06));
        assertEquals(2, negativeOffset.bucket(0.08));
        assertEquals(2, negativeOffset.bucket(0.10));
        assertEquals(3, negativeOffset.bucket(0.11));

        final WeightedPloidyHistogram positiveOffset = new WeightedPloidyHistogram(10, 0.05, 0.02);
        assertEquals(0, positiveOffset.bucket(0.04));
        assertEquals(1, positiveOffset.bucket(0.05));
        assertEquals(1, positiveOffset.bucket(0.07));
        assertEquals(1, positiveOffset.bucket(0.09));
        assertEquals(2, positiveOffset.bucket(0.10));
    }

    @Test
    public void testBucketNoOffset()
    {
        final WeightedPloidyHistogram victim = new WeightedPloidyHistogram(10, 0.01);
        assertEquals(0, victim.bucket(0.00499));
        assertEquals(1, victim.bucket(0.00500));
        assertEquals(1, victim.bucket(0.01499));
        assertEquals(2, victim.bucket(0.01500));

        assertEquals(0.00, victim.ploidy(0), 0.01);
        assertEquals(0.01, victim.ploidy(1), 0.01);
        assertEquals(0.02, victim.ploidy(2), 0.01);
    }

    @Test
    public void testHistogramConstruction()
    {
        List<WeightedPloidy> ploidies = readResource("ploidies.tsv");

        double maxPloidy = 4;
        WeightedPloidyHistogram ploidyHistogram = new WeightedPloidyHistogram(maxPloidy, 0.02);
        double[] histogram = ploidyHistogram.histogram(ploidies);

        assertEquals(3, histogram[ploidyHistogram.bucket(0.5)], 0.1);
        assertEquals(24, histogram[ploidyHistogram.bucket(1)], 0.1);
        assertEquals(8, histogram[ploidyHistogram.bucket(2)], 0.1);
        assertEquals(5, histogram[ploidyHistogram.bucket(3)], 0.1);

        ploidyHistogram = new WeightedPloidyHistogram(maxPloidy, 0.05);
        histogram = ploidyHistogram.histogram(ploidies);

        assertEquals(65, histogram[ploidyHistogram.bucket(1.0)], 0.1);
        assertEquals(20, histogram[ploidyHistogram.bucket(2.0)], 0.1);
        assertEquals(8, histogram[ploidyHistogram.bucket(3.0)], 0.1);
        assertEquals(0, histogram[ploidyHistogram.bucket(4.0)], 0.1);
    }

    @Test
    public void testHistogramPeak()
    {
        assertEquals(2, WeightedPloidyHistogram.peakBucket(1, new double[] { 10, 11, 1, 18, 8 }));
        assertEquals(3, WeightedPloidyHistogram.peakBucket(1, new double[] { 10, 11, 0, 18, 8 }));

        assertEquals(1, WeightedPloidyHistogram.peakBucket(2, new double[] { 10, 11, 0, 18, 8 }));
    }

    public static List<WeightedPloidy> readResource(@NotNull final String file)
    {
        InputStream inputStream = WeightedPloidyHistogramTest.class.getResourceAsStream("/fitting/" + file);
        List<String> lines = new BufferedReader(new InputStreamReader(inputStream)).lines().collect(Collectors.toList());

        List<WeightedPloidy> result = Lists.newArrayList();
        for(String line : lines)
        {
            String[] values = line.split(TSV_DELIM);

            WeightedPloidy ploidy = new WeightedPloidy(
                    Integer.parseInt(values[2]), Integer.parseInt(values[1]), Double.parseDouble(values[0]), 1);

            if(Doubles.lessThan(Double.parseDouble(values[0]), 10))
            {
                result.add(ploidy);
            }
        }

        return result;
    }

}
