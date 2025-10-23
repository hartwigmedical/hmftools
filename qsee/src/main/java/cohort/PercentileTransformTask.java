package cohort;

import static common.QSeeConstants.QC_LOGGER;

public class PercentileTransformTask implements Runnable
{
    private final int mFeatureIndex;

    private final FeatureMatrix mSampleFeatureMatrix;
    private PercentileTransformer mTransformer;
    private final FeatureMatrix mPercentileFeatureMatrix;

    private static final int PROGRESS_INTERVAL = 100;

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
        int numFeatures = mSampleFeatureMatrix.numFeatures();

        boolean hasFewFeatures = numFeatures <= PROGRESS_INTERVAL;
        boolean isFeatureAtInterval = mFeatureIndex % PROGRESS_INTERVAL == 0;

        if(hasFewFeatures || isFeatureAtInterval)
        {
            QC_LOGGER.debug( "Transformed {}/{} features to percentiles", mFeatureIndex, numFeatures);
        }
    }

    @Override
    public void run()
    {
        double[] featureValues = mSampleFeatureMatrix.getColumnValues(mFeatureIndex);

        mTransformer.fit(featureValues);

        String featureKey = mSampleFeatureMatrix.getFeatureKeys().get(mFeatureIndex);
        mPercentileFeatureMatrix.addColumn(featureKey, mTransformer.getRefValues());

        logProgress();

        mTransformer = null;
    }
}
