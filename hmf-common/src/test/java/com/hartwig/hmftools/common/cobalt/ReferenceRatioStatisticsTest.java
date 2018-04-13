package com.hartwig.hmftools.common.cobalt;

import static org.junit.Assert.assertEquals;

import org.jetbrains.annotations.NotNull;
import org.junit.Test;

public class ReferenceRatioStatisticsTest {

    @Test
    public void testChromosomeCount() {
        // JOBA: Klinefelter Syndrome
        assertChromosomes(true, true, 0.98, 0.48, 9842);

        // JOBA: Female
        assertChromosomes(true, false, 0.6872, 0.0049, 1304);
        assertChromosomes(true, false, 0.74, 0.01, 947);

        // JOBA: Male
        assertChromosomes(false, true, 0.50, 0.08, 9837);
        assertChromosomes(false, true, 0.51, 0.47, 9842);
        assertChromosomes(false, true, 0.51, 1.02, 9842);
    }

    private static void assertChromosomes(boolean expectedTwoX, boolean expectedY, double xMedian, double yMedian, int yCount) {
        final ReferenceRatioStatistics statistic = create(xMedian, yMedian, yCount);
        assertEquals(expectedTwoX, statistic.containsTwoXChromosomes());
        assertEquals(expectedY, statistic.containsYChromosome());
    }

    @NotNull
    private static ReferenceRatioStatistics create(double xMedian, double yMedian, int yCount) {
        return ImmutableReferenceRatioStatistics.builder().xCount(1).xMedian(xMedian).yMedian(yMedian).yCount(yCount).build();
    }
}
