package cohort;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Map;

import feature.Feature;

public class FeatureMatrix
{
    private final Map<String, double[]> mFeatureValuesMap;
    /*
    Visual representation:

             feature1     feature2     feature3  ...
    row1  array1_val1  array2_val1  array3_val1  ...
    row2  array1_val2  array2_val2  array3_val2  ...
     ...          ...          ...          ...
     */

    private final List<String> mRowIds = new ArrayList<>();
    private final int mNumRows;

    // There is no concurrent implementation of LinkedHashMap.
    // Therefore, store feature keys (= Map keys) in a list to store the insertion order of features.
    private final List<String> mFeatureKeys = new ArrayList<>();

    public FeatureMatrix(Map<String, double[]> featureValuesMap, int numRows)
    {
        if(!featureValuesMap.isEmpty())
        {
            throw new IllegalArgumentException("FeatureMatrix must be initialised with an empty map");
        }

        mFeatureValuesMap = featureValuesMap;
        mNumRows = numRows;
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
            String key = feature.mKey;

            if(!mFeatureKeys.contains(key))
            {
                mFeatureKeys.add(key);

                double[] emptyArray = new double[numRows()];
                Arrays.fill(emptyArray, Double.NaN);

                mFeatureValuesMap.put(key, emptyArray);
            }

            mFeatureValuesMap.get(key)[rowIndex] = feature.mValue;
        }
    }

    public synchronized void addColumn(String key, double[] features)
    {
        if(mFeatureKeys.contains(key))
        {
            throw new IllegalArgumentException("Cannot add a column with an already existing key");
        }

        mFeatureKeys.add(key);
        mFeatureValuesMap.put(key, features);
    }

    public void setRowIds(String[] rowIds)
    {
        if(!mRowIds.isEmpty())
        {
            throw new IllegalStateException("Cannot set row IDs if they have already been modified");
        }

        if(rowIds.length != numRows())
        {
            throw new IllegalArgumentException("Number of row IDs does not match number of initialised rows");
        }

        mRowIds.addAll(List.of(rowIds));
    }

    public int numRows() { return mNumRows; }

    public int numFeatures() { return mFeatureKeys.size(); }

    public List<String> getRowIds() { return mRowIds; }

    public List<String> getFeatureKeys() { return mFeatureKeys; }

    public double[][] getValues()
    {
        double[][] matrix = new double[numRows()][numFeatures()];

        for(int rowIndex = 0; rowIndex < numRows(); rowIndex++)
        {
            for(int featureIndex = 0; featureIndex < numFeatures(); featureIndex++)
            {
                double value = mFeatureValuesMap.get(mFeatureKeys.get(featureIndex))[rowIndex];
                matrix[rowIndex][featureIndex] = value;
            }
        }

        return matrix;
    }

    public double[][] getValuesTransposed()
    {
        double[][] matrix = new double[numFeatures()][numRows()];

        for(int rowIndex = 0; rowIndex < numRows(); rowIndex++)
        {
            for(int featureIndex = 0; featureIndex < numFeatures(); featureIndex++)
            {
                double value = mFeatureValuesMap.get(mFeatureKeys.get(featureIndex))[rowIndex];
                matrix[featureIndex][rowIndex] = value;
            }
        }

        return matrix;
    }
}
