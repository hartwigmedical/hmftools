package cohort;

import static common.QseeConstants.QC_LOGGER;

import feature.FeatureKey;

public class PercentileTransformTask implements Runnable
{
    private final int mFeatureIndex;

    private final FeatureMatrix mSampleFeatureMatrix;
    private PercentileTransformer mTransformer;
    private final FeatureMatrix mPercentileFeatureMatrix;

    public PercentileTransformTask(int featureIndex,
            FeatureMatrix sampleFeatureMatrix, PercentileTransformer transformer,
            FeatureMatrix percentileFeatureMatrix)
    {
        mFeatureIndex = featureIndex;

        mSampleFeatureMatrix = sampleFeatureMatrix;
        mTransformer = transformer;

        mPercentileFeatureMatrix = percentileFeatureMatrix;
    }

    private void logProgress()
    {
        int PROGRESS_INTERVAL = 500;
        int FEW_FEATURES_THRESHOLD = 10;

        int numFeatures = mSampleFeatureMatrix.numFeatures();

        boolean hasFewFeatures = numFeatures <= FEW_FEATURES_THRESHOLD;
        boolean isFeatureAtInterval = mFeatureIndex == numFeatures-1 || (mFeatureIndex+1) % PROGRESS_INTERVAL == 0;

        if(hasFewFeatures || isFeatureAtInterval)
        {
            QC_LOGGER.debug( "Transformed {}/{} features to percentiles", mFeatureIndex+1, numFeatures);
        }
    }

    @Override
    public void run()
    {
        double[] featureValues = mSampleFeatureMatrix.getColumnValues(mFeatureIndex);

        mTransformer.fit(featureValues);

        FeatureKey featureKey = mSampleFeatureMatrix.getFeatureKeys().get(mFeatureIndex);
        mPercentileFeatureMatrix.addColumn(featureKey, mTransformer.getRefValues());

        logProgress();

        mTransformer = null;
    }
}
