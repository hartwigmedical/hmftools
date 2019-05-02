package com.hartwig.hmftools.sig_analyser;

import static com.hartwig.hmftools.sig_analyser.common.DataUtils.greaterThan;

import static junit.framework.TestCase.assertEquals;
import static junit.framework.TestCase.assertFalse;
import static junit.framework.TestCase.assertTrue;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.sig_analyser.buckets.SigContribOptimiser;

import org.apache.logging.log4j.Level;
import org.apache.logging.log4j.core.config.Configurator;
import org.junit.Test;

public class SigContribOptimiserTests
{
    @Test
    public void testSampleFit1()
    {
        int bucketCount = 5;
        int sigCount = 3;

        List<double[]> ratiosCollection = Lists.newArrayList();

        double[] sig1 = { 0.3, 0.3, 0.05, 0.25, 0.1 };
        ratiosCollection.add(sig1);

        double[] sig2 = { 0.20, 0.10, 0.25, 0.05, 0.40 };
        ratiosCollection.add(sig2);

        double[] sig3 = { 0.0, 0.60, 0.1, 0.1, 0.2 };
        ratiosCollection.add(sig3);

        // the answer is 60, 40, 20

        double[] actualContribs = {60, 40, 20};
        double[] counts = new double[bucketCount];

        for (int j = 0; j < sigCount; ++j)
        {
            final double[] sigRatios = ratiosCollection.get(j);

            for(int i = 0; i < bucketCount; ++i)
            {
                counts[i] += actualContribs[j] * sigRatios[i];
            }
        }

        double[] noiseCounts = new double[bucketCount];

        int sample = 0;

        SigContribOptimiser sigOptim = new SigContribOptimiser(bucketCount, true, 1.0);
        sigOptim.initialise(sample, counts, noiseCounts, ratiosCollection, 0.001, 0);
        boolean calcOk = sigOptim.fitToSample();

        if (!calcOk)
            return;

        final double[] finalContribs = sigOptim.getContribs();

        // validatation that fitted counts are below the actuals + noise
        for (int i = 0; i < bucketCount; ++i)
        {
            double fittedCount = 0;

            for (int j = 0; j < sigCount; ++j)
            {
                final double[] sigRatios = ratiosCollection.get(j);
                fittedCount += sigRatios[i] * finalContribs[j];
            }

            if (greaterThan(fittedCount, counts[i] + noiseCounts[i]))
            {
                // LOGGER.error(String.format("bucket(%d) fitted count(%f) exceeds count(%f) + noise(%f)", i, fittedCount, counts[i], countsMargin[i]));
                assertFalse(true);
                break;
            }
        }

        assertEquals(sigOptim.getAllocPerc(), 1, 0.001);
    }

    @Test
    public void sampleFitTest2()
    {
        int bucketCount = 5;
        int sigCount = 4;

        double[] counts = new double[bucketCount];

        List<double[]> ratiosCollection = Lists.newArrayList();

        double[] sig1 = { 0.2, 0.2, 0.2, 0.2, 0.2 };
        ratiosCollection.add(sig1);

        double[] sig2 = { 0.0, 0.60, 0.1, 0.1, 0.20 };
        ratiosCollection.add(sig2);

        double[] sig3 = { 0.4, 0.05, 0.05, 0.4, 0.1 };
        ratiosCollection.add(sig3);

        double[] sig4 = { 0.15, 0.20, 0.25, 0.15, 0.25 };
        ratiosCollection.add(sig4);

        double[] actualContribs = { 40, 30, 60, 20 };

        for (int j = 0; j < sigCount; ++j)
        {
            final double[] sigRatios = ratiosCollection.get(j);

            for (int i = 0; i < bucketCount; ++i)
            {
                counts[i] += actualContribs[j] * sigRatios[i];
            }
        }

        // set initial contribution
        double[] countsMargin = new double[bucketCount]; // { 2, 3, 1, 1, 2 };

        // copyVector(actualContribs, contribs);
        int sample = 0;

        SigContribOptimiser sigOptim = new SigContribOptimiser(bucketCount, true, 1.0);
        sigOptim.initialise(sample, counts, countsMargin, ratiosCollection, 0.001, 0);
        boolean calcOk = sigOptim.fitToSample();

        assertTrue(calcOk);

        final double[] finalContribs = sigOptim.getContribs();

        // validatation that fitted counts are below the actuals + noise
        for (int i = 0; i < bucketCount; ++i)
        {
            double fittedCount = 0;

            for (int j = 0; j < sigCount; ++j)
            {
                final double[] sigRatios = ratiosCollection.get(j);
                fittedCount += sigRatios[i] * finalContribs[j];
            }

            if (greaterThan(fittedCount, counts[i] + countsMargin[i]))
            {
                // LOGGER.error(String.format("bucket(%d) fitted count(%f) exceeds count(%f) + noise(%f)", i, fittedCount, counts[i], countsMargin[i]));
                assertTrue(false);
                break;
            }
        }

        assertEquals(sigOptim.getAllocPerc(), 1, 0.001);
    }


    @Test
    public void testSampleFitWithNoise()
    {
        Configurator.setRootLevel(Level.DEBUG);

        int bucketCount = 5;
        int sigCount = 3;

        List<double[]> ratiosCollection = Lists.newArrayList();

        double[] sig1 = { 0.3, 0.3, 0.05, 0.25, 0.1 };
        ratiosCollection.add(sig1);

        double[] sig2 = { 0.20, 0.10, 0.25, 0.05, 0.40 };
        ratiosCollection.add(sig2);

        double[] sig3 = { 0.0, 0.60, 0.1, 0.1, 0.2 };
        ratiosCollection.add(sig3);

        // the answer is 60, 40, 20

        double[] actualContribs = {60, 40, 20};
        double[] counts = new double[bucketCount];

        for (int j = 0; j < sigCount; ++j)
        {
            final double[] sigRatios = ratiosCollection.get(j);

            for(int i = 0; i < bucketCount; ++i)
            {
                counts[i] += actualContribs[j] * sigRatios[i];
            }
        }

        double[] noiseCounts = new double[bucketCount];
        double bucketNoise = 5;
        for(int i = 0; i < bucketCount; ++i)
        {
            noiseCounts[i] = bucketNoise;
        }

        int sample = 0;

        SigContribOptimiser sigOptim = new SigContribOptimiser(bucketCount, true, 1.0);
        sigOptim.initialise(sample, counts, noiseCounts, ratiosCollection, 0.001, 0);
        boolean calcOk = sigOptim.fitToSample();

        if (!calcOk)
            return;

        final double[] finalContribs = sigOptim.getContribs();

        // validatation that fitted counts are below the actual counts + noise
        for (int i = 0; i < bucketCount; ++i)
        {
            double fittedCount = 0;

            for (int j = 0; j < sigCount; ++j)
            {
                final double[] sigRatios = ratiosCollection.get(j);
                fittedCount += sigRatios[i] * finalContribs[j];
            }

            if (greaterThan(fittedCount, counts[i] + noiseCounts[i]))
            {
                // LOGGER.error(String.format("bucket(%d) fitted count(%f) exceeds count(%f) + noise(%f)", i, fittedCount, counts[i], countsMargin[i]));
                assertTrue(false);
                break;
            }
        }

        assertEquals(sigOptim.getAllocPerc(), 1, 0.001);
    }

}
