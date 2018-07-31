package com.hartwig.hmftools.data_analyser;

import org.apache.commons.math3.distribution.ChiSquaredDistribution;
import org.apache.commons.math3.distribution.PoissonDistribution;

public class MiscTester {


    public static void runTests()
    {
        chiSquaredTests();

    }

    private static void chiSquaredTests()
    {
        int degreesOfFreedom = 95;
        ChiSquaredDistribution chiSquDist = new ChiSquaredDistribution(degreesOfFreedom);

        double result = chiSquDist.cumulativeProbability(135);
    }

    private static void poissonTests()
    {
        PoissonDistribution poisson = new PoissonDistribution(100);

        double prob = poisson.cumulativeProbability(99);
        prob = poisson.cumulativeProbability(110);
    }

}
