package com.hartwig.hmftools.purple.fitting;

import static org.junit.Assert.assertEquals;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.purple.fittingsnv.SomaticKernelDensityPeaks;
import com.hartwig.hmftools.purple.fittingsnv.SomaticPeak;

import org.jetbrains.annotations.NotNull;
import org.junit.Test;

public class SomaticKernelDensityPeaksTest
{

    private static final double EPSILON = 1e-10;

    @Test
    public void testFindPeaks()
    {
        List<Double> sample = Lists.newArrayList();
        for(int i = 0; i < 100; i++)
        {
            sample.add(i / 100d);
        }

        for(int i = 0; i < 10; i++)
        {
            sample.add(0.2);
            sample.add(0.4);
            sample.add(0.6);
        }

        final List<SomaticPeak> peaks = SomaticKernelDensityPeaks.findPeaks(sample);
        assertEquals(2, peaks.size());
        assertPeak(peaks.get(0), 0.2, 13);
        assertPeak(peaks.get(1), 0.4, 13);
    }

    private static void assertPeak(@NotNull SomaticPeak victim, double vaf, int count)
    {
        assertEquals(vaf, victim.AlleleFrequency, EPSILON);
        assertEquals(count, victim.Count);
    }
}


