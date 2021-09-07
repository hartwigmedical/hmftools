package com.hartwig.hmftools.geneutils.cooc;

import com.hartwig.hmftools.common.stats.FisherExactTest;

import org.junit.Test;

public class FisherTest
{
    @Test
    public void testStatsRoutines()
    {
        FisherExactTest fetCalc = new FisherExactTest();
        fetCalc.initialise(1000);

        int withAwithB = 11;
        int withANoB = 27;
        int noAWithB = 2;
        int noAnoB = 170;
        double expectedCount = 5;

        double fisherProb = fetCalc.getLeftTailedP(withAwithB, noAWithB, withANoB, noAnoB);
        fisherProb = fetCalc.getRightTailedP(withAwithB, noAWithB, withANoB, noAnoB);

    }
}
