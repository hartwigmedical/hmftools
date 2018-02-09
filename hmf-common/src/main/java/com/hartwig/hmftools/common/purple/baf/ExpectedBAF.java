package com.hartwig.hmftools.common.purple.baf;

import org.apache.commons.math3.distribution.BinomialDistribution;

public class ExpectedBAF {

    public static double expectedBAF(int averageDepth) {

        int minDepth = averageDepth / 2;
        int maxDepth = averageDepth * 3 / 2;
        double totalProbability = 0;

        int j = 0;
        for (int i = minDepth; i < maxDepth; i++) {
            final BinomialDistribution distribution = new BinomialDistribution(i, 0.5);
            double probability = 1d * distribution.inverseCumulativeProbability(0.75) / i;
            totalProbability += probability;
            j++;
        }

        return totalProbability / (maxDepth - minDepth);
    }

}
