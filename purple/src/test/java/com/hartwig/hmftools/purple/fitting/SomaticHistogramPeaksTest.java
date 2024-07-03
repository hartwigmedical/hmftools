package com.hartwig.hmftools.purple.fitting;

import static java.lang.Math.round;

import static com.google.common.collect.Lists.newArrayList;

import static junit.framework.TestCase.assertEquals;

import java.util.List;

import com.hartwig.hmftools.purple.fittingsnv.SomaticHistogramPeaks;
import com.hartwig.hmftools.purple.fittingsnv.WeightedPloidy;

import org.junit.Test;

public class SomaticHistogramPeaksTest
{
    @Test
    public void testFindPeaks()
    {
        final List<WeightedPloidy> weightedVAFs = newArrayList();

        addVafEntry(weightedVAFs, 0.04, 1, 3);
        addVafEntry(weightedVAFs, 0.05, 1, 6);
        addVafEntry(weightedVAFs, 0.06, 1, 8);
        addVafEntry(weightedVAFs, 0.07, 1, 11);
        addVafEntry(weightedVAFs, 0.08, 1, 7);
        addVafEntry(weightedVAFs, 0.09, 1, 3);
        addVafEntry(weightedVAFs, 0.10, 1, 5);
        addVafEntry(weightedVAFs, 0.11, 1, 3);
        addVafEntry(weightedVAFs, 0.12, 1, 4);
        addVafEntry(weightedVAFs, 0.13, 1, 2);
        addVafEntry(weightedVAFs, 0.14, 1, 4);
        addVafEntry(weightedVAFs, 0.15, 1, 2);
        addVafEntry(weightedVAFs, 0.16, 1, 1);

        SomaticHistogramPeaks somaticHistogram = new SomaticHistogramPeaks(2, 0.01, 1, 12, 0.1);
        double maxPurity = somaticHistogram.findVafPeak(weightedVAFs);
        assertEquals(0.07, maxPurity, 0.001);

        // relax the conditions to select a higher peak
        somaticHistogram = new SomaticHistogramPeaks(2, 0.01, 1, 5, 0.05);
        maxPurity = somaticHistogram.findVafPeak(weightedVAFs);
        assertEquals(0.15, maxPurity, 0.001);

        // threshold too high so use the max-weighted peak
        somaticHistogram = new SomaticHistogramPeaks(2, 0.01, 1, 50, 0.05);
        maxPurity = somaticHistogram.findVafPeak(weightedVAFs);
        assertEquals(0.07, maxPurity, 0.001);
    }

    private static void addVafEntry(final List<WeightedPloidy> weightedPloidies, double vaf, double weight, int multiples)
    {
        int totalReadCount = 100;
        int alleleReadCount = (int) round(vaf * totalReadCount);

        for(int i = 0; i < multiples; ++i)
        {
            weightedPloidies.add(new WeightedPloidy(totalReadCount, alleleReadCount, vaf, weight));
        }
    }
}
