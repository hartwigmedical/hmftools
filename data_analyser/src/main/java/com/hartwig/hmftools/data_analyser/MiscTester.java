package com.hartwig.hmftools.data_analyser;

import static com.hartwig.hmftools.data_analyser.calcs.DataUtils.copyVector;
import static com.hartwig.hmftools.data_analyser.calcs.NmfSampleFitter.fitCountsToRatios;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.data_analyser.calcs.NmfSampleFitter;
import com.hartwig.hmftools.data_analyser.calcs.SigContributionOptimiser;

import org.apache.commons.math3.distribution.ChiSquaredDistribution;
import org.apache.commons.math3.distribution.PoissonDistribution;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

public class MiscTester {

    private static final Logger LOGGER = LogManager.getLogger(MiscTester.class);

    public static void runTests()
    {
        sampleFitTest();

        // chiSquaredTests();

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

    private static void sampleFitTest()
    {
        int bucketCount = 5;
        int sigCount = 3;

        // double[] counts = new double[bucketCount];

        List<double[]> ratiosCollection = Lists.newArrayList();

        double[] sig1 = { 0.3, 0.3, 0.05, 0.25, 0.1 };
        ratiosCollection.add(sig1);

        double[] sig2 = { 0.20, 0.10, 0.25, 0.05, 0.40 };
        ratiosCollection.add(sig2);

        double[] sig3 = { 0.0, 0.60, 0.1, 0.1, 0.2 };
        ratiosCollection.add(sig3);

        // the answer is 60, 40, 20

        /*
        double[] actualContribs = {60, 40, 20};

        for (int j = 0; j < sigCount; ++j)
        {
            final double[] sigRatios = ratiosCollection.get(j);

            for(int i = 0; i < bucketCount; ++i)
            {
                counts[i] += actualContribs[j] * sigRatios[i];
            }
        }
        */

        double[] counts = {26, 34, 15, 19, 26};

//        counts[0] = 35;
//        counts[1] = 32;
//        counts[2] = 16;
//        counts[3] = 18;
//        counts[4] = 24;

        // set initial contribution
        double[] contribs = new double[sigCount];

        double[] countsMargin = new double[bucketCount]; // { 2, 3, 1, 1, 2 };

        // copyVector(actualContribs, contribs);

        contribs[0] = 70;
        contribs[1] = 50;
        contribs[2] = 30;

        int sample = 0;

        // boolean calcOk = fitCountsToRatios(sample, counts, countsMargin, ratiosCollection, contribs, 0.001);

        SigContributionOptimiser sigOptim = new SigContributionOptimiser(sample, counts, countsMargin, ratiosCollection, contribs);
        boolean calcOk = sigOptim.fitToSample(0.001);

        // boolean calcOk = fitCountsToRatios_v2(sample, counts, countsMargin, ratiosCollection, contribs, 0.001);

        if (!calcOk)
            return;

        final double[] finalContribs = sigOptim.getContribs();

        // validatation that fitted counts are below the actuals + noise
        boolean allOk = true;
        for (int i = 0; i < bucketCount; ++i)
        {
            double fittedCount = 0;

            for (int j = 0; j < sigCount; ++j)
            {
                final double[] sigRatios = ratiosCollection.get(j);
                fittedCount += sigRatios[i] * finalContribs[j];
            }

            if (fittedCount > counts[i] + countsMargin[i])
            {
                LOGGER.error(String.format("bucket(%d) fitted count(%.1f) exceeds count(%.0f) + noise(%.0f)", i, fittedCount, counts[i], countsMargin[i]));
                allOk = false;
            }
        }

        if(allOk)
            LOGGER.debug("sig optimisation test success");
    }

}
