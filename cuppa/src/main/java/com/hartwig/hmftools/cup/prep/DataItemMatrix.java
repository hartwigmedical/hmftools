package com.hartwig.hmftools.cup.prep;

import static com.hartwig.hmftools.cup.common.CupConstants.CUP_LOGGER;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import java.util.Map;

import com.google.common.annotations.VisibleForTesting;

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

        checkDimensions();
    }

    private void checkDimensions()
    {
        for(String[] row : FeatureBySampleMatrix.values())
        {
            if(nSamples() != row.length)
            {
                CUP_LOGGER.error("Found FeatureBySampleMatrix row with length {}, but required {} (= number of SampleIds)", row.length, nSamples());
                System.exit(1);
            }
        }
    }

    public String[] get(DataItem.Index index)
    {
        return FeatureBySampleMatrix.get(index);
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
        Collections.sort(Indexes);
    }

    public void printRows()
    {
        for(DataItem.Index index : Indexes)
        {
            System.out.println(index + ", Values=" + Arrays.toString(get(index)));
        }
    }

    @VisibleForTesting
    public String[] getSampleFeatureValues(int sampleIndex)
    {
        String[] featureValues = new String[nFeatures()];

        int i = 0;
        for(DataItem.Index dataItemIndex : Indexes)
        {
            featureValues[i] = get(dataItemIndex)[sampleIndex];
            i++;
        }

        return featureValues;
    }
}
