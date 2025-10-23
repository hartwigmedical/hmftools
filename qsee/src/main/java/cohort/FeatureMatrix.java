package cohort;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;

import feature.Feature;
import feature.FeatureKey;

public class FeatureMatrix
{
    private final Map<FeatureKey, double[]> mFeatureValuesMap;
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
    private final List<FeatureKey> mFeatureKeys = new ArrayList<>();

    private static final double EMPTY_VALUE = Double.NaN;

    public FeatureMatrix(Map<FeatureKey, double[]> featureValuesMap, int numRows)
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
            FeatureKey key = feature.key();
            addColumnIfMissing(key);
            mFeatureValuesMap.get(key)[rowIndex] = feature.value();
        }
    }

    private void addColumnIfMissing(FeatureKey key)
    {
        if(!mFeatureKeys.contains(key))
        {
            mFeatureKeys.add(key);

            double[] emptyArray = new double[numRows()];
            Arrays.fill(emptyArray, EMPTY_VALUE);

            mFeatureValuesMap.put(key, emptyArray);
        }
    }

    public synchronized void addColumn(FeatureKey key, double[] features)
    {
        if(mFeatureKeys.contains(key))
        {
            throw new IllegalArgumentException("Cannot add a column with an already existing key");
        }

        mFeatureKeys.add(key);
        mFeatureValuesMap.put(key, features);
    }

    public void setRowIds(List<String> rowIds)
    {
        if(!mRowIds.isEmpty())
        {
            throw new IllegalStateException("Cannot set row IDs if they have already been modified");
        }

        if(rowIds.size() != numRows())
        {
            throw new IllegalArgumentException("No. of row IDs does not match number of initialised rows");
        }

        mRowIds.addAll(rowIds);
    }

    public FeatureMatrix reorderRows(List<String> rowIdsOrdered)
    {
        if(rowIdsOrdered.size() != numRows())
        {
            throw new IllegalArgumentException(String.format("Mismatched no. of rows (%d) provided and no. of rows in matrix (%d)",
                    rowIdsOrdered.size(), numRows()
            ));
        }

        if(!new HashSet<>(rowIdsOrdered).containsAll(mRowIds))
        {
            throw new IllegalArgumentException("Row IDs provided do not match the row IDs in the matrix");
        }

        FeatureMatrix newMatrix = new FeatureMatrix(new HashMap<>(), numRows());
        newMatrix.setRowIds(rowIdsOrdered);

        for(int newRowIndex = 0; newRowIndex < numRows(); newRowIndex++)
        {
            String newRowId = rowIdsOrdered.get(newRowIndex);
            int oldRowIndex = mRowIds.indexOf(newRowId);

            for(int featureIndex = 0; featureIndex < numFeatures(); featureIndex++)
            {
                FeatureKey featureKey = mFeatureKeys.get(featureIndex);
                double value = mFeatureValuesMap.get(featureKey)[oldRowIndex];

                newMatrix.addColumnIfMissing(featureKey);
                newMatrix.mFeatureValuesMap.get(featureKey)[newRowIndex] = value;
            }
        }

        return newMatrix;
    }

    public int numRows() { return mNumRows; }

    public int numFeatures() { return mFeatureKeys.size(); }

    public List<String> getRowIds() { return mRowIds; }

    public List<FeatureKey> getFeatureKeys() { return mFeatureKeys; }

    public double[][] getValues()
    {
        double[][] matrix = new double[numRows()][numFeatures()];

        for(int rowIndex = 0; rowIndex < numRows(); rowIndex++)
        {
            matrix[rowIndex] = getRowValues(rowIndex);
        }

        return matrix;
    }

    public double[][] getValuesTransposed()
    {
        double[][] matrix = new double[numFeatures()][numRows()];

        for(int featureIndex = 0; featureIndex < numFeatures(); featureIndex++)
        {
            matrix[featureIndex] = getColumnValues(mFeatureKeys.get(featureIndex));
        }

        return matrix;
    }

    public double[] getRowValues(int index)
    {
        double[] rowValues = new double[numFeatures()];

        for(int featureIndex = 0; featureIndex < numFeatures(); featureIndex++)
        {
            rowValues[featureIndex] = mFeatureValuesMap.get(mFeatureKeys.get(featureIndex))[index];
        }

        return rowValues;
    }

    public double[] getRowValues(String rowId)
    {
        int rowIndex = mRowIds.indexOf(rowId);

        if(rowIndex < 0)
        {
            throw new IllegalArgumentException(String.format("RowId(%s) not found", rowId));
        }

        return getRowValues(rowIndex);
    }

    public double[] getColumnValues(int index)
    {
        return getColumnValues(mFeatureKeys.get(index));
    }

    public double[] getColumnValues(FeatureKey key)
    {
        return mFeatureValuesMap.get(key);
    }
}
