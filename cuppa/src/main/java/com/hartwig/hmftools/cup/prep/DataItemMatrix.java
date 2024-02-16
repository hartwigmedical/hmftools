package com.hartwig.hmftools.cup.prep;

public class DataItemMatrix
{
    public final String[] SampleIds;
    public final DataSource[] Sources;
    public final ItemType[] Types;
    public final String[] Keys;
    public final String[][] FeatureBySampleMatrix;

    public DataItemMatrix(
            final String[] sampleIds,
            final DataSource[] sources,
            final ItemType[] types,
            final String[] keys,
            final String[][] featureBySampleMatrix
    ){
        SampleIds = sampleIds;

        Sources = sources;
        Types = types;
        Keys = keys;

        FeatureBySampleMatrix = featureBySampleMatrix;

        checkRows();
    }

    public int getNFeatures()
    {
        return FeatureBySampleMatrix.length;
    }

    public int getNSamples()
    {
        return FeatureBySampleMatrix[0].length;
    }

    private void checkRows()
    {
        int expectedLength = getNFeatures();

        int[] arrayLengths = {
                Sources.length,
                Types.length,
                Keys.length,
                FeatureBySampleMatrix.length
        };

        for(int arrayLength : arrayLengths)
        {
            if(arrayLength != expectedLength)
                throw new IllegalStateException("All input arrays must have the same size");
        }
    }
}
