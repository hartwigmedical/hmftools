package com.hartwig.hmftools.sigs;

import org.junit.Test;

public class SigSampleTest
{
    @Test
    public void testSampleAllocActuals()
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
}
