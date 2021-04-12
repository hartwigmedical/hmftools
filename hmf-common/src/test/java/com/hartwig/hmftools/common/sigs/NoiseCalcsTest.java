package com.hartwig.hmftools.common.sigs;

import static com.hartwig.hmftools.common.sigs.NoiseCalcs.calcPoissonRangeGivenProb;

import static org.junit.Assert.assertEquals;

import org.junit.Test;

public class NoiseCalcsTest
{
    @Test
    public void testPoissonRange()
    {
        double probability = 0.0001;

        int range = calcPoissonRangeGivenProb(10, probability, 10, true);
        assertEquals(10, range);

        range = calcPoissonRangeGivenProb(20, probability, 10, true);
        assertEquals(15, range);

        range = calcPoissonRangeGivenProb(100, probability, 10, true);
        assertEquals(36, range);

        range = calcPoissonRangeGivenProb(1000, probability, 10, true);
        assertEquals(116, range);

        range = calcPoissonRangeGivenProb(10000, probability, 10, true);
        assertEquals(370, range);

        range = calcPoissonRangeGivenProb(100000, probability, 10, true);
        assertEquals(1170, range);

        range = calcPoissonRangeGivenProb(1000000, probability, 10, true);
        assertEquals(3700, range);

        // now test finding the upper value
        range = calcPoissonRangeGivenProb(1, probability, 0, false);
        assertEquals(5, range);

        range = calcPoissonRangeGivenProb(2, probability, 0, false);
        assertEquals(7, range);

        range = calcPoissonRangeGivenProb(5, probability, 0, false);
        assertEquals(10, range);

        range = calcPoissonRangeGivenProb(10, probability, 0, false);
        assertEquals(14, range);

        range = calcPoissonRangeGivenProb(100, probability, 0, false);
        assertEquals(40, range);
    }

    /*
    // chiSquaredTests
    int degreesOfFreedom = 95;
    ChiSquaredDistribution chiSquDist = new ChiSquaredDistribution(degreesOfFreedom);

    double result = chiSquDist.cumulativeProbability(135);

    PoissonDistribution poisson = new PoissonDistribution(100);

    double prob = poisson.cumulativeProbability(99);
    prob = poisson.cumulativeProbability(110);

    */


}
