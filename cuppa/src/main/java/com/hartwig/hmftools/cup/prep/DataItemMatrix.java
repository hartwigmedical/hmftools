package com.hartwig.hmftools.cup.prep;

import java.util.concurrent.ConcurrentHashMap;

public class DataItemMatrix
{
    public final String[] SampleIds;
    public final ConcurrentHashMap<DataItem.Index, String[]> FeatureBySampleMatrix;

    public DataItemMatrix(
            final String[] sampleIds,
            final ConcurrentHashMap<DataItem.Index, String[]> featureBySampleMatrix
    ){
        SampleIds = sampleIds;
        FeatureBySampleMatrix = featureBySampleMatrix;
    }

    public String[] get(DataItem.Index index)
    {
        return FeatureBySampleMatrix.get(index);
    }

    public void put(DataItem.Index index, String[] values)
    {
        FeatureBySampleMatrix.put(index, values);
    }

    public DataItem.Index[] getFeatureIndexes()
    {
        return FeatureBySampleMatrix.keySet().toArray(DataItem.Index[]::new);
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
