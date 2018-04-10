package com.hartwig.hmftools.common.cobalt;

import static org.junit.Assert.assertEquals;

import org.junit.Test;

public class ReferenceRatioStatisticsTest {
    @Test
    public void testChromosomeCount() {

        // Klinefelter Syndrome
        assertChromsomes(true, true, 0.98, 0.48, 9842);

        // Female
        assertChromsomes(true, false, 0.6872, 0.0049, 1304);
        assertChromsomes(true, false, 0.74, 0.01, 947);

        // Male
        assertChromsomes(false, true, 0.50, 0.08, 9837);
        assertChromsomes(false, true, 0.51, 0.47, 9842);
        assertChromsomes(false, true, 0.51, 1.02, 9842);
    }

    private void assertChromsomes(boolean expectedTwoX, boolean expectedY,  double xMedian, double yMedian, int yCount) {
        final ReferenceRatioStatistics statistic = create(xMedian, yMedian, yCount);
        assertEquals(expectedTwoX, statistic.containsTwoXChromosomes());
        assertEquals(expectedY, statistic.containsYChromosome());
    }


    private ReferenceRatioStatistics create(double xMedian, double yMedian, int yCount) {
        return ImmutableReferenceRatioStatistics.builder().xCount(1).xMedian(xMedian).yMedian(yMedian).yCount(yCount).build();
    }

}
