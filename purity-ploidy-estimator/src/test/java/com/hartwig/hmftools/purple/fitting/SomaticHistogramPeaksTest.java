package com.hartwig.hmftools.purple.fitting;

import static java.lang.Math.round;

import static com.google.common.collect.Lists.newArrayList;

import static junit.framework.TestCase.assertEquals;

import java.util.List;

import com.hartwig.hmftools.common.purple.purity.SomaticPeak;
import com.hartwig.hmftools.common.variant.clonality.ModifiableWeightedPloidy;

import org.apache.commons.math3.distribution.BinomialDistribution;
import org.junit.Test;

public class SomaticHistogramPeaksTest
{
    @Test
    public void testFindPeaks()
    {
        final List<ModifiableWeightedPloidy> weightedVAFs = newArrayList();

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

        double maxPurity = SomaticHistogramPeaks.findVafPeak(2, 0.01, 1, weightedVAFs);
        assertEquals(0.07, maxPurity, 0.001);
    }

    private static void addVafEntry(final List<ModifiableWeightedPloidy> weightedPloidies, double vaf, double weight, int multiples)
    {
        int totalReadCount = 100;
        int alleleReadCount = (int)round(vaf * totalReadCount);

        for(int i = 0; i < multiples; ++i)
        {
            weightedPloidies.add(ModifiableWeightedPloidy.create()
                    .setPloidy(vaf)
                    .setWeight(weight)
                    .setAlleleReadCount(alleleReadCount)
                    .setTotalReadCount(totalReadCount));
        }
    }
}
