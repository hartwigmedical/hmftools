package cohort;

import java.util.ArrayList;
import java.util.List;
import java.util.Map;

import feature.FeatureValue;

public class FeatureMatrix
{
    private final Map<String, Double[]> mFeatureValuesMap;

    private final List<String> mRowIds;
    private final int mNumRows;

    // There is no concurrent implementation of LinkedHashMap.
    // Therefore, store feature keys (= Map keys) in a list to store the insertion order of features.
    private final List<String> mFeatureKeys;

    public FeatureMatrix(Map<String, Double[]> featureValuesMap, int numRows)
    {
        if(!featureValuesMap.isEmpty())
        {
            throw new IllegalArgumentException("SampleFeatureMatrix must be initialised with an empty map");
        }

        mFeatureValuesMap = featureValuesMap;

        mRowIds = new ArrayList<>();
        mNumRows = numRows;

        mFeatureKeys = new ArrayList<>();
    }

    public synchronized void addFeaturesToRow(String rowId, List<FeatureValue> features)
    {
        if(!mRowIds.contains(rowId))
        {
            mRowIds.add(rowId);
        }

        int rowIndex = mRowIds.indexOf(rowId);

        for(FeatureValue feature : features)
        {
            String key = feature.mKey;

            if(!mFeatureKeys.contains(key))
            {
                mFeatureKeys.add(key);
                mFeatureValuesMap.put(key, new Double[numRows()]);
            }

            mFeatureValuesMap.get(key)[rowIndex] = feature.mValue;
        }
    }

    public int numRows() { return mNumRows; }

    public int numFeatures() { return mFeatureKeys.size(); }

    public List<String> getRowIds() { return mRowIds; }

    public List<String> getFeatureKeys() { return mFeatureKeys; }

    public double[][] getValues(double nullFillValue)
    {
        double[][] matrix = new double[numRows()][numFeatures()];

        for(int rowIndex = 0; rowIndex < numRows(); rowIndex++)
        {
            for(int featureIndex = 0; featureIndex < numFeatures(); featureIndex++)
            {
                Double sampleFeatureValue = mFeatureValuesMap.get(mFeatureKeys.get(featureIndex))[rowIndex];
                matrix[rowIndex][featureIndex] = sampleFeatureValue != null ? sampleFeatureValue : nullFillValue;
            }
        }

        return matrix;
    }

    public double[][] getValuesTransposed(double nullFillValue)
    {
        double[][] matrix = new double[numFeatures()][numRows()];

        for(int rowIndex = 0; rowIndex < numRows(); rowIndex++)
        {
            for(int featureIndex = 0; featureIndex < numFeatures(); featureIndex++)
            {
                Double sampleFeatureValue = mFeatureValuesMap.get(mFeatureKeys.get(featureIndex))[rowIndex];
                matrix[featureIndex][rowIndex] = sampleFeatureValue != null ? sampleFeatureValue : nullFillValue;
            }
        }

        return matrix;
    }
}
