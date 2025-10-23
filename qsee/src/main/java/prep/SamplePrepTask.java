package prep;

import static common.QseeConstants.QC_LOGGER;

import java.util.ArrayList;
import java.util.List;

import org.jetbrains.annotations.Nullable;

import cohort.FeatureMatrix;
import feature.Feature;

public class SamplePrepTask implements Runnable
{
    private final PrepConfig mConfig;
    private final int mSampleIndex;
    private final List<CategoryPrep> mCategoryPreps;

    private final List<Feature> mFeatures = new ArrayList<>();

    private static final int PROGRESS_INTERVAL = 100;

    @Nullable
    private final FeatureMatrix mSampleFeatureMatrix;

    public SamplePrepTask(PrepConfig config, int sampleIndex, List<CategoryPrep> categoryPreps,
            @Nullable FeatureMatrix sampleFeatureMatrix)
    {
        mConfig = config;
        mSampleIndex = sampleIndex;
        mCategoryPreps = categoryPreps;

        mSampleFeatureMatrix = sampleFeatureMatrix;
    }

    private void logProgress(int sampleIndex, int categoryIndex)
    {
        int sampleCount = mConfig.SampleIds.size();
        String sampleId = mConfig.SampleIds.get(sampleIndex);

        int categoryCount = mCategoryPreps.size();
        String categoryName = mCategoryPreps.get(categoryIndex).getClass().getSimpleName();

        if(sampleCount == 1)
        {
            QC_LOGGER.debug("Extracting data for sample({}}) task({}/{}: {})",
                    sampleId, categoryIndex, categoryCount, categoryName);
        }
        else
        {
            boolean hasFewSamples = sampleCount <= PROGRESS_INTERVAL;
            boolean isSampleAtInterval = sampleIndex == 0 || sampleIndex == sampleCount-1 ||
                    (sampleIndex+1) % PROGRESS_INTERVAL == 0;

            if(hasFewSamples || isSampleAtInterval)
            {
                QC_LOGGER.debug("Extracting data for sample({}/{}: {}) task({}/{}: {})",
                        sampleIndex+1, sampleCount, sampleId,
                        categoryIndex+1, categoryCount, categoryName);
            }
        }
    }

    @Override
    public void run()
    {
        String sampleId = mConfig.SampleIds.get(mSampleIndex);

        for(int categoryIndex = 0; categoryIndex < mCategoryPreps.size(); categoryIndex++)
        {
            logProgress(mSampleIndex, categoryIndex);

            CategoryPrep categoryPrep = mCategoryPreps.get(categoryIndex);
            List<Feature> categoryFeatures = categoryPrep.extractSampleData(sampleId);
            mFeatures.addAll(categoryFeatures);
        }

        if(mSampleFeatureMatrix != null)
        {
            mSampleFeatureMatrix.addRow(sampleId, mFeatures);
        }

        mFeatures.clear();
    }
}