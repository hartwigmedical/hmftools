package com.hartwig.hmftools.sigs.buckets;

import static com.hartwig.hmftools.common.sigs.DataUtils.greaterThan;
import static com.hartwig.hmftools.common.sigs.DataUtils.lessOrEqual;

import static junit.framework.TestCase.assertEquals;
import static junit.framework.TestCase.assertTrue;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.sigs.LeastSquaresFit;
import com.hartwig.hmftools.common.utils.Matrix;
import com.hartwig.hmftools.sigs.fitter.ConstrainedFitter;

import org.apache.logging.log4j.Level;
import org.apache.logging.log4j.core.config.Configurator;
import org.junit.Ignore;
import org.junit.Test;

public class SigContribOptimiserTest
{
    @Ignore
    @Test
    public void testSampleFit1()
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

        {
            Matrix sigs = new Matrix(bucketCount, sigCount);
            sigs.setCol(0, sig1);
            sigs.setCol(1, sig2);
            sigs.setCol(2, sig3);

            LeastSquaresFit lsqFit = new LeastSquaresFit(bucketCount, sigCount);
            lsqFit.initialise(sigs.getData(), counts);

            lsqFit.solve();
            final double[] sigAllocs = lsqFit.getContribs();
            assertEquals(40, sigAllocs[0], 0.001);
        }


        ConstrainedFitter fitter = new ConstrainedFitter(bucketCount);
        fitter.setParameters(0, null, true, 1.0);
        boolean calcOk = fitter.fitCounts(ratiosCollection, counts, null);

        assertTrue(calcOk);

        final double[] finalContribs = fitter.getContribs();

        // validation that fitted counts are below the actuals + noise
        for (int i = 0; i < bucketCount; ++i)
        {
            double fittedCount = 0;

            for (int j = 0; j < sigCount; ++j)
            {
                final double[] sigRatios = ratiosCollection.get(j);
                fittedCount += sigRatios[i] * finalContribs[j];
            }

            assertTrue(lessOrEqual(fittedCount, counts[i]));
        }

        // assertEquals(1, fitter.getAllocPerc(), 0.01);
    }

    @Ignore
    @Test
    public void testSampleFit_Old1()
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

        int sampleId = 0;
        SampleData sample = new SampleData(sampleId);
        sample.setBucketCounts(counts);
        sample.setElevatedBucketCounts(counts, noiseCounts);

        SampleSigContribOptimiser sigOptim = new SampleSigContribOptimiser(bucketCount, true, 1.0);
        sigOptim.initialise(sample, ratiosCollection, 0.001, 0);
        sigOptim.setMinContribChange(0.01, 0.001);
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

            assertTrue(lessOrEqual(fittedCount, counts[i] + noiseCounts[i]));
        }

        assertEquals(1, sigOptim.getAllocPerc(), 0.01);
    }

    @Ignore
    @Test
    public void testSampleFit2()
    {
        Configurator.setRootLevel(Level.DEBUG);

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

        {
            Matrix sigs = new Matrix(bucketCount, sigCount);
            sigs.setCol(0, sig1);
            sigs.setCol(1, sig2);
            sigs.setCol(2, sig3);
            sigs.setCol(3, sig4);

            LeastSquaresFit lsqFit = new LeastSquaresFit(bucketCount, sigCount);
            lsqFit.initialise(sigs.getData(), counts);

            lsqFit.solve();
            final double[] sigAllocs = lsqFit.getContribs();
            assertEquals(40, sigAllocs[0], 0.001);
        }

        double[] noiseCounts = new double[bucketCount];

        int sampleId = 0;
        SampleData sample = new SampleData(sampleId);
        sample.setBucketCounts(counts);
        sample.setElevatedBucketCounts(counts, noiseCounts);

        SampleSigContribOptimiser sigOptim = new SampleSigContribOptimiser(bucketCount, true, 1.0);
        sigOptim.initialise(sample, ratiosCollection, 0.001, 0);
        boolean calcOk = sigOptim.fitToSample();

        assertTrue(calcOk);

        final double[] finalContribs = sigOptim.getContribs();

        // validation that fitted counts are below the actuals + noise
        for (int i = 0; i < bucketCount; ++i)
        {
            double fittedCount = 0;

            for (int j = 0; j < sigCount; ++j)
            {
                final double[] sigRatios = ratiosCollection.get(j);
                fittedCount += sigRatios[i] * finalContribs[j];
            }

            assertTrue(lessOrEqual(fittedCount, counts[i] + noiseCounts[i]));
        }

        assertEquals(1, sigOptim.getAllocPerc(), 0.003);
    }

    @Test
    @Ignore
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

        int sampleId = 0;
        SampleData sample = new SampleData(sampleId);
        sample.setBucketCounts(counts);
        sample.setElevatedBucketCounts(counts, noiseCounts);

        SampleSigContribOptimiser sigOptim = new SampleSigContribOptimiser(bucketCount, true, 1.0);
        // sigOptim.initialise(sampleId, counts, noiseCounts, ratiosCollection, 0.001, 0);
        sigOptim.initialise(sample, ratiosCollection, 0.001, 0);
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

        assertEquals(1, sigOptim.getAllocPerc(), 0.001);
    }

}
