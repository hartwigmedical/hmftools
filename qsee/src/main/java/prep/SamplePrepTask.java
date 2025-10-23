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
    private final CategoryPrep mCategoryPrep;

    private final List<Feature> mFeatures = new ArrayList<>();

    @Nullable
    private final FeatureMatrix mSampleFeatureMatrix;

    public SamplePrepTask(PrepConfig config, int sampleIndex, CategoryPrep categoryPrep,
            @Nullable FeatureMatrix sampleFeatureMatrix)
    {
        mConfig = config;
        mSampleIndex = sampleIndex;
        mCategoryPrep = categoryPrep;

        mSampleFeatureMatrix = sampleFeatureMatrix;
    }

    private void logProgress(int sampleIndex)
    {
        int PROGRESS_INTERVAL = 100;
        int FEW_SAMPLES_THRESHOLD = 10;

        int sampleCount = mConfig.SampleIds.size();
        String sampleId = mConfig.SampleIds.get(sampleIndex);

        if(sampleCount == 1)
        {
            QC_LOGGER.debug("Extracting data for sample: {}}", sampleId);
        }
        else
        {
            boolean hasFewSamples = sampleCount <= FEW_SAMPLES_THRESHOLD;
            boolean isSampleAtInterval = sampleIndex == 0 || sampleIndex == sampleCount-1 ||
                    (sampleIndex+1) % PROGRESS_INTERVAL == 0;

            if(hasFewSamples || isSampleAtInterval)
            {
                QC_LOGGER.debug("Extracting data for sample {}/{}: {}", sampleIndex+1, sampleCount, sampleId);
            }
        }
    }

    @Override
    public void run()
    {
        String sampleId = mConfig.SampleIds.get(mSampleIndex);

        logProgress(mSampleIndex);

        List<Feature> categoryFeatures = mCategoryPrep.extractSampleData(sampleId);
        mFeatures.addAll(categoryFeatures);

        if(mSampleFeatureMatrix != null)
        {
            mSampleFeatureMatrix.addRow(sampleId, mFeatures);
            mFeatures.clear();
        }
    }
}