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
        mFeatureValuesMap = featureValuesMap;

        mSampleIds = new ArrayList<>();
        mNumSamples = numSamples;

        mFeatureKeys = new ArrayList<>();
    }

//    public SampleFeatureMatrix(Map<String, Double[]> featureValuesMap, List<String> sampleIds)
//    {
//        mFeatureValuesMap = featureValuesMap;
//
//        mSampleIds = sampleIds;
//        mNumSamples = sampleIds.size();
//
//        mFeatureKeys = new ArrayList<>(featureValuesMap.keySet());
//
//        validateInput();
//    }
//
//    private void validateInput()
//    {
//        if(mSampleIds.isEmpty())
//        {
//            throw new IllegalArgumentException("Sample IDs must not be empty");
//        }
//
//        Set<String> uniqueSampleIds = new HashSet<>(mSampleIds);
//        int nonUniqueSampleIdCount = mSampleIds.size() - uniqueSampleIds.size();
//
//        if(nonUniqueSampleIdCount > 0)
//        {
//            throw new IllegalArgumentException(String.format("Sample IDs must be unique. Found %s non-unique sample IDs",
//                    nonUniqueSampleIdCount));
//        }
//
//        if(!mFeatureValuesMap.isEmpty())
//        {
//            String firstKey = mFeatureValuesMap.keySet().stream().findFirst().get();
//            int expectedNumSamples = mFeatureValuesMap.get(firstKey).length;
//
//            for(String key : mFeatureValuesMap.keySet())
//            {
//                if(mFeatureValuesMap.get(key).length != expectedNumSamples)
//                {
//                    throw new IllegalArgumentException(String.format("All feature columns must have the same length. Expected %s but got %s for key '%s'",
//                            expectedNumSamples, mFeatureValuesMap.get(key).length, key));
//                }
//            }
//        }
//    }

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
                mFeatureValuesMap.put(key, new Double[getNumSamples()]);
            }

            mFeatureValuesMap.get(key)[sampleIndex] = feature.mValue;
        }
    }

    public int getNumSamples() { return mNumSamples; }

    public int getNumFeatures() { return mFeatureKeys.size(); }

    public List<String> getSampleIds(){ return mSampleIds; }

    public List<String> getFeatureKeys(){ return mFeatureKeys; }

    public Double[][] getValues(boolean samplesAsRows)
    {
        return samplesAsRows ? getSampleByFeature2DArray() : getFeatureBySample2DArray();
    }

    private Double[][] getSampleByFeature2DArray()
    {
        Double[][] matrix = new Double[getNumSamples()][getNumFeatures()];

        for(int i = 0; i < getNumSamples(); i++)
        {
            for(int j = 0; j < getNumSamples(); j++)
            {
                matrix[i][j] = mFeatureValuesMap.get(mFeatureKeys.get(j))[i];
            }
        }

        return matrix;
    }

    private Double[][] getFeatureBySample2DArray()
    {
        Double[][] matrix = new Double[getNumFeatures()][getNumSamples()];

        for(int i = 0; i < getNumFeatures(); i++)
        {
            matrix[i] = mFeatureValuesMap.get(mFeatureKeys.get(i));
        }

        return matrix;
    }
}
