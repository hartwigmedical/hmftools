package com.hartwig.hmftools.data_analyser;

import org.apache.commons.math3.distribution.PoissonDistribution;

public class MiscTester {


    public static void runTests()
    {

        PoissonDistribution poisson = new PoissonDistribution(100);

        double prob = poisson.cumulativeProbability(99);
        prob = poisson.cumulativeProbability(110);

    }
}
