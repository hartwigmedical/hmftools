package com.hartwig.hmftools.cup.prep;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import java.util.Map;

public class DataItemMatrix
{
    public final String[] SampleIds;
    public final Map<DataItem.Index, String[]> FeatureBySampleMatrix;
    public List<DataItem.Index> Indexes;

    public DataItemMatrix(
            final String[] sampleIds,
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

    public void put(DataItem.Index index, String[] values)
    {
        FeatureBySampleMatrix.put(index, values);
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
        return SampleIds.length;
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
}
