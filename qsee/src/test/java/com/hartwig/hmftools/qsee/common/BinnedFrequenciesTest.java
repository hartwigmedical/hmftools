package com.hartwig.hmftools.qsee.common;

import static org.junit.Assert.assertArrayEquals;

import org.junit.Test;

public class BinnedFrequenciesTest
{
    private static final long[] BIN_STARTS    = { 0  , 10 , 25 , 50 , 100 };
    private static final double[] FREQUENCIES = { 200, 200, 200, 200, 200 };

    private static final BinnedFrequencies BINNED_FREQUENCIES = new BinnedFrequencies(BIN_STARTS, FREQUENCIES);

    @Test
    public void canCalcCorrectBinEnds()
    {
        long[] binEnds = BINNED_FREQUENCIES.binEnds();
        long[] binEndsExpected = { 10, 25, 50, 100, Long.MAX_VALUE };
        assertArrayEquals(binEndsExpected, binEnds);
    }

    @Test
    public void canCalcFrequencyDensities()
    {
        double[] frequencyDensities = BINNED_FREQUENCIES.calcFrequencyDensities();
        double[] frequencyDensitiesExpected = {
                (double) 200 / 10,
                (double) 200 / (25 - 10),
                (double) 200 / (50 - 25),
                (double) 200 / (100 - 50),
                200
        };
        assertArrayEquals(frequencyDensitiesExpected, frequencyDensities, 0.01);
    }

    @Test
    public void canCalcProportionalDensities(){
        double[] proportionalDensities = BINNED_FREQUENCIES.calcProportionalDensities();
        double[] proportionalDensitiesExpected = { 0.02, 0.01334, 0.008, 0.004, 0.2 };
        assertArrayEquals(proportionalDensitiesExpected, proportionalDensities, 0.001);
    }
}
