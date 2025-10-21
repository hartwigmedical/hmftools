package cohort;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import feature.FeatureValue;

public class SampleFeatureMatrix
{
    private final Map<String, Double[]> mFeatureValuesMap;
    private final List<String> mSampleIds;

    // There is no concurrent implementation of LinkedHashMap.
    // Therefore, store feature keys (= Map keys) in a list to store the insertion order of features.
    private final List<String> mFeatureKeys = new ArrayList<>();

    public SampleFeatureMatrix(Map<String, Double[]> featureValuesMap, List<String> sampleIds)
    {
        mFeatureValuesMap = featureValuesMap;
        mSampleIds = sampleIds;

        validateInput();

        mFeatureKeys.addAll(featureValuesMap.keySet());
    }

    private void validateInput()
    {
        if(mSampleIds.isEmpty())
        {
            throw new IllegalArgumentException("Sample IDs must not be empty");
        }

        Set<String> uniqueSampleIds = new HashSet<>(mSampleIds);
        int nonUniqueSampleIdCount = mSampleIds.size() - uniqueSampleIds.size();

        if(nonUniqueSampleIdCount > 0)
        {
            throw new IllegalArgumentException(String.format("Sample IDs must be unique. Found %s non-unique sample IDs",
                    nonUniqueSampleIdCount));
        }

        if(!mFeatureValuesMap.isEmpty())
        {
            String firstKey = mFeatureValuesMap.keySet().stream().findFirst().get();
            int expectedNumSamples = mFeatureValuesMap.get(firstKey).length;

            for(String key : mFeatureValuesMap.keySet())
            {
                if(mFeatureValuesMap.get(key).length != expectedNumSamples)
                {
                    throw new IllegalArgumentException(String.format("All feature columns must have the same length. Expected %s but got %s for key '%s'",
                            expectedNumSamples, mFeatureValuesMap.get(key).length, key));
                }
            }
        }
    }

    public List<String> getSampleIds(){ return mSampleIds; }

    public List<String> getFeatureKeys(){ return mFeatureKeys; }

    public Double[][] getValues(boolean samplesAsRows)
    {
        return samplesAsRows ? getSampleByFeatureMatrix() : getFeatureBySampleMatrix();
    }

    private Double[][] getSampleByFeatureMatrix()
    {
        Double[][] matrix = new Double[mSampleIds.size()][mFeatureKeys.size()];

        for(int i = 0; i < mSampleIds.size(); i++)
        {
            for(int j = 0; j < mFeatureKeys.size(); j++)
            {
                matrix[i][j] = mFeatureValuesMap.get(mFeatureKeys.get(j))[i];
            }
        }

        return matrix;
    }

    private Double[][] getFeatureBySampleMatrix()
    {
        Double[][] matrix = new Double[mFeatureKeys.size()][mSampleIds.size()];

        for(int i = 0; i < mFeatureKeys.size(); i++)
        {
            matrix[i] = mFeatureValuesMap.get(mFeatureKeys.get(i));
        }

        return matrix;
    }

    public synchronized void addSampleFeatures(String sampleId, List<FeatureValue<Double>> features)
    {
        if(!mSampleIds.contains(sampleId))
        {
            mSampleIds.add(sampleId);
        }

        int sampleIndex = mSampleIds.indexOf(sampleId);

        for(FeatureValue<Double> feature : features)
        {
            String key = feature.mKey;

            if(!mFeatureKeys.contains(key))
            {
                mFeatureKeys.add(key);
                mFeatureValuesMap.put(key, new Double[mSampleIds.size()]);
            }

            mFeatureValuesMap.get(key)[sampleIndex] = feature.mValue;
        }
    }
}
