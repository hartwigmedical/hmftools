package com.hartwig.hmftools.data_analyser;

import static com.hartwig.hmftools.data_analyser.calcs.DataUtils.copyVector;
import static com.hartwig.hmftools.data_analyser.calcs.DataUtils.sumVector;
import static com.hartwig.hmftools.data_analyser.calcs.NmfSampleFitter.fitCountsToRatios;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.data_analyser.calcs.DataUtils;
import com.hartwig.hmftools.data_analyser.calcs.SigContributionOptimiser;
import com.hartwig.hmftools.data_analyser.loaders.GenericDataLoader;
import com.hartwig.hmftools.data_analyser.types.GenericDataCollection;
import com.hartwig.hmftools.data_analyser.types.NmfMatrix;
import com.hartwig.hmftools.data_analyser.types.SampleData;

import org.apache.commons.math3.distribution.ChiSquaredDistribution;
import org.apache.commons.math3.distribution.PoissonDistribution;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

public class MiscTester
{

    private static final Logger LOGGER = LogManager.getLogger(MiscTester.class);

    public static void runTests()
    {
        //sampleFitTest();
        //sampleFitTest2();
        testSampleAllocActuals();
        // testSigOptimiserActuals();

        // stringTest();
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
        LOGGER.info("sampleFitTest");

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

        double[] counts = { 26, 34, 15, 19, 26 };

        //        counts[0] = 35;
        //        counts[1] = 32;
        //        counts[2] = 16;
        //        counts[3] = 18;
        //        counts[4] = 24;

        // set initial contribution
        double[] contribs = new double[sigCount];

        double[] countsMargin = new double[bucketCount]; // { 2, 3, 1, 1, 2 };

        // copyVector(actualContribs, contribs);

        int sample = 0;

        // boolean calcOk = fitCountsToRatios(sample, counts, countsMargin, ratiosCollection, contribs, 0.001);

        SigContributionOptimiser sigOptim = new SigContributionOptimiser(bucketCount, true, 1.0, false);
        sigOptim.initialise(sample, counts, countsMargin, ratiosCollection, 0.001, 0);
        boolean calcOk = sigOptim.fitToSample();

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

        if (!allOk)
            LOGGER.debug("sampleFitTest: error");

        if (sigOptim.getAllocPerc() == 1)
            LOGGER.debug("sampleFitTest: success");
        else
            LOGGER.debug("sampleFitTest: non-optimal fit");
    }

    private static void sampleFitTest2()
    {
        LOGGER.info("sampleFitTest2");

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

        // double[] counts = {26, 34, 15, 19, 26};

        // set initial contribution
        double[] countsMargin = new double[bucketCount]; // { 2, 3, 1, 1, 2 };

        // copyVector(actualContribs, contribs);
        int sample = 0;

        // boolean calcOk = fitCountsToRatios(sample, counts, countsMargin, ratiosCollection, contribs, 0.001);

        SigContributionOptimiser sigOptim = new SigContributionOptimiser(bucketCount, true, 1.0, false);
        sigOptim.initialise(sample, counts, countsMargin, ratiosCollection, 0.001, 0);
        boolean calcOk = sigOptim.fitToSample();

        if (!calcOk)
        {
            LOGGER.error("sampleFitTest2 failed");
            return;
        }

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

        if (!allOk)
            LOGGER.debug("sampleFitTest2: error");

        if (sigOptim.getAllocPerc() == 1)
            LOGGER.debug("sampleFitTest2: success");
        else
            LOGGER.debug("sampleFitTest2: non-optimal fit");
    }

    private static void stringTest()
    {
        String tmp = "AID=0.23131;UV=0.93";

        String tmp2 = tmp.replaceAll("[0-9.=]", "");
    }

    private static void testSampleAllocActuals()
    {
        // load counts and sigs to fit with
        String countsFile = "/Users/charlesshale/data/r_data/snv_nmf_matrix_data.csv";
        GenericDataCollection dataCollection = GenericDataLoader.loadFile(countsFile);
        NmfMatrix sampleCountsMatrix = DataUtils.createMatrixFromListData(dataCollection.getData());
        sampleCountsMatrix.cacheTranspose();

        String noiseFile = "/Users/charlesshale/dev/nmf/snv_ba_sample_noise.csv";
        dataCollection = GenericDataLoader.loadFile(noiseFile);
        NmfMatrix sampleNoiseMatrix = DataUtils.createMatrixFromListData(dataCollection.getData());
        sampleNoiseMatrix.cacheTranspose();

        String sigsFile = "/Users/charlesshale/dev/nmf/snv_ba_predefined_sigs.csv";
        dataCollection = GenericDataLoader.loadFile(sigsFile);
        NmfMatrix sigs = DataUtils.createMatrixFromListData(dataCollection.getData());
        sigs.cacheTranspose();

        int sampleId = 2207;

        final double[] sampleCounts = sampleCountsMatrix.getCol(sampleId);
        final double[] sampleNoise = sampleNoiseMatrix.getCol(sampleId);

        SampleData sample = new SampleData(sampleId);
        sample.setBucketCounts(sampleCounts);
        sample.setElevatedBucketCounts(sampleCounts, sampleNoise);

        // try adding the specific groups in turn
        double[] bgRatios = sigs.getCol(6);

        List<Integer> bucketIds = Lists.newArrayList();
        for(int b = 0; b < bgRatios.length; ++b)
        {
            if(bgRatios[b] > 0)
                bucketIds.add(b);
        }

        double[] sampleAllocCounts = sample.getPotentialElevCounts(bgRatios, bucketIds);
        double allocTotal = sumVector(sampleAllocCounts);
        sample.allocateBucketCounts(sampleAllocCounts, 0.03);

        bgRatios = sigs.getCol(22);

        bucketIds.clear();
        for(int b = 0; b < bgRatios.length; ++b)
        {
            if(bgRatios[b] > 0)
                bucketIds.add(b);
        }

        // sampleAllocCounts = sample.getPotentialElevCounts(bgRatios, bucketIds);
        sampleAllocCounts = sample.getPotentialUnallocCounts(bgRatios, bucketIds);
        allocTotal = sumVector(sampleAllocCounts);
        sample.allocateBucketCounts(sampleAllocCounts, 0.03);

    }

    private static void testSigOptimiserActuals()
    {
        // load counts and sigs to fit with
        String countsFile = "/Users/charlesshale/data/r_data/snv_nmf_matrix_data.csv";
        GenericDataCollection dataCollection = GenericDataLoader.loadFile(countsFile);
        NmfMatrix sampleCountsMatrix = DataUtils.createMatrixFromListData(dataCollection.getData());
        sampleCountsMatrix.cacheTranspose();

        String noiseFile = "/Users/charlesshale/dev/nmf/snv_ba_sample_noise.csv";
        dataCollection = GenericDataLoader.loadFile(noiseFile);
        NmfMatrix sampleNoiseMatrix = DataUtils.createMatrixFromListData(dataCollection.getData());
        sampleNoiseMatrix.cacheTranspose();

        String sigsFile = "/Users/charlesshale/dev/nmf/snv_ba_predefined_sigs.csv";
        dataCollection = GenericDataLoader.loadFile(sigsFile);
        NmfMatrix sigs = DataUtils.createMatrixFromListData(dataCollection.getData());
        sigs.cacheTranspose();


        int sampleId = 2207;

        final double[] sampleCounts = sampleCountsMatrix.getCol(sampleId);
        final double[] sampleNoise = sampleNoiseMatrix.getCol(sampleId);

        SampleData sample = new SampleData(sampleId);
        sample.setBucketCounts(sampleCounts);
        sample.setElevatedBucketCounts(sampleCounts, sampleNoise);

        List<double[]> ratiosCollection = Lists.newArrayList();

        ratiosCollection.add(sigs.getCol(22)); // bg
        ratiosCollection.add(sigs.getCol(6)); // bg
        /*
        ratiosCollection.add(sigs.getCol(42)); // bg
        ratiosCollection.add(sigs.getCol(43));
        ratiosCollection.add(sigs.getCol(12));
        ratiosCollection.add(sigs.getCol(37));
        */

        double[] contribs = new double[ratiosCollection.size()];

        SigContributionOptimiser sigOptim = new SigContributionOptimiser(sampleCounts.length, true, 0.999, true);
        sigOptim.initialise(sampleId, sampleCounts, sampleNoise, ratiosCollection, 0.03, 400);
        sigOptim.setRequiredSig(1);

        List<Integer> sigIds = Lists.newArrayList();
        sigIds.add(22);
        sigIds.add(6);
        /*
        sigIds.add(1471);
        sigIds.add(1502);
        sigIds.add(12);
        sigIds.add(1450);
        */
        sigOptim.setSigIds(sigIds);

        boolean calcOk = sigOptim.fitToSample();

        if (!calcOk)
            return;

        // final double[] finalContribs = sigOptim.getContribs();


    }

}
