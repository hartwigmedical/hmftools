package cohort;

import java.util.ArrayList;
import java.util.List;
import java.util.Map;

import feature.FeatureValue;

public class SampleFeatureMatrix
{
    private final Map<String, Double[]> mFeatureValuesMap;

    private final List<String> mSampleIds;
    private final int mNumSamples;

    // There is no concurrent implementation of LinkedHashMap.
    // Therefore, store feature keys (= Map keys) in a list to store the insertion order of features.
    private final List<String> mFeatureKeys;

    public SampleFeatureMatrix(Map<String, Double[]> featureValuesMap, int numSamples)
    {
        if(!featureValuesMap.isEmpty())
        {
            throw new IllegalArgumentException("SampleFeatureMatrix must be initialised with an empty map");
        }

        mFeatureValuesMap = featureValuesMap;

        mSampleIds = new ArrayList<>();
        mNumSamples = numSamples;

        mFeatureKeys = new ArrayList<>();
    }

    public synchronized void addSampleFeatures(String sampleId, List<FeatureValue> features)
    {
        if(!mSampleIds.contains(sampleId))
        {
            mSampleIds.add(sampleId);
        }

        int sampleIndex = mSampleIds.indexOf(sampleId);

        for(FeatureValue feature : features)
        {
            String key = feature.mKey;

            if(!mFeatureKeys.contains(key))
            {
                mFeatureKeys.add(key);
                mFeatureValuesMap.put(key, new Double[numSamples()]);
            }

            mFeatureValuesMap.get(key)[sampleIndex] = feature.mValue;
        }
    }

    public int numSamples() { return mNumSamples; }

    public int numFeatures() { return mFeatureKeys.size(); }

    public List<String> getSampleIds() { return mSampleIds; }

    public List<String> getFeatureKeys() { return mFeatureKeys; }

    public double[][] getValues(double nullFillValue)
    {
        double[][] matrix = new double[numSamples()][numFeatures()];

        for(int sampleIndex = 0; sampleIndex < numSamples(); sampleIndex++)
        {
            for(int featureIndex = 0; featureIndex < numFeatures(); featureIndex++)
            {
                Double sampleFeatureValue = mFeatureValuesMap.get(mFeatureKeys.get(featureIndex))[sampleIndex];
                matrix[sampleIndex][featureIndex] = sampleFeatureValue != null ? sampleFeatureValue : nullFillValue;
            }
        }

        return matrix;
    }
}
