package com.hartwig.hmftools.cup.prep;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import java.util.Map;

public class DataItemMatrix
{
    public final List<String> SampleIds;
    public final Map<DataItem.Index, String[]> FeatureBySampleMatrix;
    public List<DataItem.Index> Indexes;

    public DataItemMatrix(
            final List<String> sampleIds,
            final Map<DataItem.Index, String[]> featureBySampleMatrix
    ){
        SampleIds = sampleIds;
        FeatureBySampleMatrix = featureBySampleMatrix;
        Indexes = new ArrayList<>(FeatureBySampleMatrix.keySet());
    }

    public String[] get(DataItem.Index index)
    {
        return FeatureBySampleMatrix.get(index);
    }

    public List<DataItem.Index> getIndexes()
    {
        return Indexes;
    }

    public int nFeatures()
    {
        return Indexes.size();
    }

    public int nSamples()
    {
        return SampleIds.size();
    }

    public void sortIndexes()
    {
        Collections.sort(Indexes, new DataItem.IndexComparator());
    }

    public void printRows()
    {
        for(DataItem.Index index : Indexes)
        {
            System.out.println(index + " Values=" + Arrays.toString(get(index)));
        }
    }

    public List<String> getFeatureValuesBySampleIndex(int sampleIndex)
    {
        List<String> featureValues = new ArrayList<>();
        for(DataItem.Index featureIndex : Indexes)
        {
            featureValues.add(get(featureIndex)[sampleIndex]);
        }
        return featureValues;
    }
}
