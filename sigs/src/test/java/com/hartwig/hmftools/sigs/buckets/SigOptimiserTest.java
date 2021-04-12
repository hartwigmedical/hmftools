package com.hartwig.hmftools.sigs.buckets;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.junit.Test;

public class SigOptimiserTest
{
    private static final Logger LOGGER = LogManager.getLogger(SigOptimiserTest.class);


    @Test
    public void testCalcRatiosAndRanges()
    {
        /*
        // load counts and sigs to fit with
        String countsFile = "/Users/charlesshale/data/r_data/snv_nmf_matrix_data.csv";
        GenericDataCollection dataCollection = GenericDataLoader.loadFile(countsFile);
        SigMatrix sampleCountsMatrix = DataUtils.createMatrixFromListData(dataCollection.getData());
        sampleCountsMatrix.cacheTranspose();

        String noiseFile = "/Users/charlesshale/dev/nmf/snv_ba_sample_noise.csv";
        dataCollection = GenericDataLoader.loadFile(noiseFile);
        SigMatrix sampleNoiseMatrix = DataUtils.createMatrixFromListData(dataCollection.getData());
        sampleNoiseMatrix.cacheTranspose();

        String sigsFile = "/Users/charlesshale/dev/nmf/snv_ba_predefined_sigs.csv";
        dataCollection = GenericDataLoader.loadFile(sigsFile);
        SigMatrix sigs = DataUtils.createMatrixFromListData(dataCollection.getData());
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

        double[] sampleAllocCounts = sample.getPotentialCounts(bgRatios, bucketIds, null);
        double allocTotal = sumVector(sampleAllocCounts);
        sample.allocateBucketCounts(sampleAllocCounts, 0.03);

        bgRatios = sigs.getCol(22);

        bucketIds.clear();
        for(int b = 0; b < bgRatios.length; ++b)
        {
            if(bgRatios[b] > 0)
                bucketIds.add(b);
        }

        // sampleAllocCounts = sample.getPotentialCounts(bgRatios, bucketIds);
        sampleAllocCounts = sample.getPotentialUnallocCounts(bgRatios, bucketIds, null);
        allocTotal = sumVector(sampleAllocCounts);
        sample.allocateBucketCounts(sampleAllocCounts, 0.03);
        */

    }

    @Test
    public void testSigOptimiserActuals()
    {
        /*
        // load counts and sigs to fit with
        String countsFile = "/Users/charlesshale/data/r_data/snv_nmf_matrix_data.csv";
        GenericDataCollection dataCollection = GenericDataLoader.loadFile(countsFile);
        SigMatrix sampleCountsMatrix = DataUtils.createMatrixFromListData(dataCollection.getData());
        sampleCountsMatrix.cacheTranspose();

        String noiseFile = "/Users/charlesshale/dev/nmf/snv_ba_sample_noise.csv";
        dataCollection = GenericDataLoader.loadFile(noiseFile);
        SigMatrix sampleNoiseMatrix = DataUtils.createMatrixFromListData(dataCollection.getData());
        sampleNoiseMatrix.cacheTranspose();

        String sigsFile = "/Users/charlesshale/dev/nmf/snv_ba_predefined_sigs.csv";
        dataCollection = GenericDataLoader.loadFile(sigsFile);
        SigMatrix sigs = DataUtils.createMatrixFromListData(dataCollection.getData());
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

        double[] contribs = new double[ratiosCollection.size()];

        SigContribOptimiser sigOptim = new SigContribOptimiser(sampleCounts.length, true, 0.999);
        sigOptim.initialise(sampleId, sampleCounts, sampleNoise, ratiosCollection, 0.03, 400);
        sigOptim.setRequiredSig(1);

        List<Integer> sigIds = Lists.newArrayList();
        sigIds.add(22);
        sigIds.add(6);
        sigOptim.setSigIds(sigIds);

        boolean calcOk = sigOptim.fitToSample();

        if (!calcOk)
            return;

        // final double[] finalContribs = sigOptim.getContribs();
        */

    }

    @Test
    public void testSigReconstruction()
    {
        /*
        String sigsFile = "/Users/charlesshale/dev/nmf/snv_ba_predefined_sigs.csv";
        GenericDataCollection dataCollection = GenericDataLoader.loadFile(sigsFile);
        dataCollection = GenericDataLoader.loadFile(sigsFile);
        SigMatrix sigs = DataUtils.createMatrixFromListData(dataCollection.getData());
        sigs.cacheTranspose();

        int bgSigCount = 20;

        int bucketCount = sigs.Rows;

        SigContribOptimiser sigOptim = new SigContribOptimiser(bucketCount, false, 0.99);

        double[] testGroupRatios = new double[bucketCount];
        double[] testGroupNoise = new double[bucketCount]; // unused
        List<double[]> ratiosCollection = Lists.newArrayList();
        List<Integer> sigIds = Lists.newArrayList();
        List<Integer> testSigBuckets = Lists.newArrayList();
        double minRatioThreshold = 0.001;

        for(int testSigId = bgSigCount; testSigId < sigs.Cols; ++testSigId)
        {
            copyVector(sigs.getCol(testSigId), testGroupRatios);
            vectorMultiply(testGroupRatios, 10000); // to make the other group contributions be a percentage of total
            testSigBuckets.clear();

            for(int b = 0; b < bucketCount; ++b)
            {
                if(testGroupRatios[b] > 0)
                    testSigBuckets.add(b);
            }

            ratiosCollection.clear();
            sigIds.clear();

            for(int otherSigId = bgSigCount; otherSigId < sigs.Cols; ++otherSigId)
            {
                if (otherSigId == testSigId)
                    continue;

                double[] otherGroupRatios = new double[bucketCount];
                copyVector(sigs.getCol(otherSigId), otherGroupRatios);

                boolean hasDiffBuckets = false;
                for(int b = 0; b < bucketCount; ++b)
                {
                    if(otherGroupRatios[b] > 0 && !testSigBuckets.contains(b))
                    {
                        if(otherGroupRatios[b] < minRatioThreshold)
                        {
                            otherGroupRatios[b] = 0; // skip this bucket and continue on
                            continue;
                        }

                        hasDiffBuckets = true;
                        break;
                    }
                }

                if(hasDiffBuckets)
                    continue;

                // vectorMultiply(otherGroupRatios, 100);

                ratiosCollection.add(otherGroupRatios);
                sigIds.add(otherSigId);
            }

            if(ratiosCollection.size() < 2)
            {
                LOGGER.debug(String.format("bg(%d) insufficient overlapping sigs for reconstruction", testSigId));
                continue;
            }

            sigOptim.initialise(testSigId, testGroupRatios, testGroupNoise, ratiosCollection, 0.01, 0);
            sigOptim.setSigIds(sigIds);
            sigOptim.setLogVerbose(true);

            boolean validCalc = sigOptim.fitToSample();

            if(!validCalc)
                continue;

            if(sigOptim.getAllocPerc() < 0.95)
            {
                LOGGER.debug(String.format("bg(%d) achieved low reconstruction from %d sigs to %.3f percent",
                        testSigId, sigOptim.contributingSigCount(), sigOptim.getAllocPerc()));
                continue;
            }

            LOGGER.debug(String.format("bg(%d) achieved high reconstruction from %d sigs to %.3f percent:",
                    testSigId, sigOptim.contributingSigCount(), sigOptim.getAllocPerc()));

            final double[] sigContribs = sigOptim.getContribs();
            for(int sig = 0; sig < ratiosCollection.size(); ++sig)
            {
                LOGGER.debug(String.format("bg(%d) from sigig(%d) contrib(%.3f) percent", testSigId, sig, sigContribs[sig]/100));
            }
        }
        */
    }

}
