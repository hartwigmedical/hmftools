package com.hartwig.hmftools.qsee.cohort;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.stream.Stream;

import com.google.common.annotations.VisibleForTesting;
import com.hartwig.hmftools.qsee.feature.Feature;
import com.hartwig.hmftools.qsee.feature.FeatureKey;

public class FeatureMatrix
{
    private final Map<FeatureKey, Feature[]> mFeatureValuesMap;
    /*
    Visual representation:

             feature1  feature2  feature3  ...
    sample1     s1_f1     s1_f2     s1_f3  ...
    sample2     s2_f1     s2_f2     s2_f3  ...
        ...       ...       ...       ...
     */

    private final List<String> mRowIds = new ArrayList<>();
    private final int mNumRows;

    // There is no concurrent implementation of LinkedHashMap.
    // Therefore, store feature keys (= Map keys) in a list to store the insertion order of features.
    private final List<FeatureKey> mFeatureKeys = new ArrayList<>();

    public FeatureMatrix(Map<FeatureKey, Feature[]> featureValuesMap, List<String> expectedRowIds)
    {
        mFeatureValuesMap = featureValuesMap;
        checkMapEmpty();

        mRowIds.addAll(expectedRowIds);
        mNumRows = expectedRowIds.size();
    }

    public FeatureMatrix(Map<FeatureKey, Feature[]> featureValuesMap, int numRows)
    {
        mFeatureValuesMap = featureValuesMap;
        checkMapEmpty();

        mNumRows = numRows;
    }

    private void checkMapEmpty()
    {
        if(!mFeatureValuesMap.isEmpty())
        {
            throw new IllegalStateException("FeatureMatrix must be initialised with an empty map");
        }
    }

    public synchronized void addRow(String rowId, List<Feature> features)
    {
        if(!mRowIds.contains(rowId))
        {
            mRowIds.add(rowId);
        }

        int rowIndex = mRowIds.indexOf(rowId);

        for(Feature feature : features)
        {
            FeatureKey key = feature.key();
            addColumnIfMissing(key);
            mFeatureValuesMap.get(key)[rowIndex] = feature;
        }
    }

    private void addColumnIfMissing(FeatureKey key)
    {
        if(!mFeatureKeys.contains(key))
        {
            mFeatureKeys.add(key);

            Feature emptyFeature = new Feature(key, Double.NaN);
            Feature[] emptyArray = new Feature[numRows()];
            Arrays.fill(emptyArray, emptyFeature);

            mFeatureValuesMap.put(key, emptyArray);
        }
    }

    public void sortFeatureKeys()
    {
        Collections.sort(mFeatureKeys);
    }

    public int numRows() { return mNumRows; }

    public int numFeatures() { return mFeatureKeys.size(); }

    public List<String> getRowIds() { return mRowIds; }

    public List<FeatureKey> getFeatureKeys() { return mFeatureKeys; }

    public Feature[] getRow(int index)
    {
        Feature[] row = new Feature[numFeatures()];

        for(int featureIndex = 0; featureIndex < numFeatures(); featureIndex++)
        {
            row[featureIndex] = mFeatureValuesMap.get(mFeatureKeys.get(featureIndex))[index];
        }

        return row;
    }

    public Feature[] getRow(String rowId)
    {
        int rowIndex = mRowIds.indexOf(rowId);

        if(rowIndex < 0)
        {
            throw new IllegalArgumentException(String.format("RowId(%s) not found", rowId));
        }

        return getRow(rowIndex);
    }

    public Feature[] getColumn(int index)
    {
        return getColumn(mFeatureKeys.get(index));
    }

    public Feature[] getColumn(FeatureKey key)
    {
        return mFeatureValuesMap.get(key);
    }

    @VisibleForTesting
    double[][] getFeatureValues()
    {
        double[][] matrix = new double[numRows()][numFeatures()];

        for(int rowIndex = 0; rowIndex < numRows(); rowIndex++)
        {
            Feature[] features = getRow(rowIndex);
            double[] featureValues = Stream.of(features).mapToDouble(Feature::value).toArray();
            matrix[rowIndex] = featureValues;
        }

        return matrix;
    }
}
