package prep;

import static common.QSeeConstants.QC_LOGGER;

import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.Callable;

import feature.FeatureValue;

public class SamplePrepTask<T> implements Callable<List<FeatureValue<T>>>
{
    private final PrepConfig mConfig;
    private final int mSampleIndex;
    private final List<CategoryPrep<T>> mCategoryPreps;

    public SamplePrepTask(PrepConfig config, int sampleIndex, List<CategoryPrep<T>> categoryPreps)
    {
        mConfig = config;
        mCategoryPreps = categoryPreps;
        mSampleIndex = sampleIndex;
    }

    @Override
    public List<FeatureValue<T>> call()
    {
        String sampleId = mConfig.SampleIds.get(mSampleIndex);

        List<FeatureValue<T>> featureValues = new ArrayList<>();

        for(CategoryPrep<T> categoryPrep : mCategoryPreps)
        {
            QC_LOGGER.debug("Running {} for sample ({}/{}): {}",
                    categoryPrep.getClass().getSimpleName(),
                    mSampleIndex + 1, mConfig.SampleIds.size(), sampleId);

            List<FeatureValue<T>> categoryFeatureValues = categoryPrep.extractSampleData(sampleId);

            featureValues.addAll(categoryFeatureValues);
        }

        return featureValues;
    }
}