package com.hartwig.hmftools.cup.prep;

import java.util.concurrent.ConcurrentHashMap;

public class DataItemMatrix
{
    public final String[] SampleIds;
    public final ConcurrentHashMap<DataItemIndex, String[]> FeatureBySampleMatrix;

    public DataItemMatrix(
            final String[] sampleIds,
            final ConcurrentHashMap<DataItemIndex, String[]> featureBySampleMatrix
    ){
        SampleIds = sampleIds;
        FeatureBySampleMatrix = featureBySampleMatrix;
    }

    public String[] get(DataItemIndex index)
    {
        return FeatureBySampleMatrix.get(index);
    }

    public void put(DataItemIndex index, String[] values)
    {
        FeatureBySampleMatrix.put(index, values);
    }

    public DataItemIndex[] getFeatureIndexes()
    {
        return FeatureBySampleMatrix.keySet().toArray(DataItemIndex[]::new);
    }

    public int nFeatures()
    {
        return FeatureBySampleMatrix.size();
    }

    public int nSamples()
    {
        return SampleIds.length;
    }
}
