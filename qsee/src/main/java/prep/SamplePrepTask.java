package prep;

import static common.QseeConstants.QC_LOGGER;

import java.util.ArrayList;
import java.util.List;

import org.jetbrains.annotations.Nullable;

import cohort.FeatureMatrix;
import feature.Feature;

public class SamplePrepTask implements Runnable
{
    private final List<String> mSampleIds;
    private final int mSampleIndex;
    private final CategoryPrep mCategoryPrep;

    private List<Feature> mFeatures;

    @Nullable
    private final FeatureMatrix mSampleFeatureMatrix;

    public SamplePrepTask(List<String> sampleIds, int sampleIndex, CategoryPrep categoryPrep,
            @Nullable FeatureMatrix sampleFeatureMatrix)
    {
        mSampleIds = sampleIds;
        mSampleIndex = sampleIndex;
        mCategoryPrep = categoryPrep;

        mSampleFeatureMatrix = sampleFeatureMatrix;
    }

    private void logProgress(int sampleIndex)
    {
        int PROGRESS_INTERVAL = 100;
        int FEW_SAMPLES_THRESHOLD = 10;

        int sampleCount = mSampleIds.size();
        String sampleId = mSampleIds.get(sampleIndex);

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
        String sampleId = mSampleIds.get(mSampleIndex);

        try
        {
            logProgress(mSampleIndex);

            mFeatures = mCategoryPrep.extractSampleData(sampleId);

            if(mSampleFeatureMatrix != null)
            {
                mSampleFeatureMatrix.addRow(sampleId, mFeatures);
                mFeatures = null;
            }
        }
        catch(Exception e)
        {
            QC_LOGGER.error("Failed to run {} for sample: {}", mCategoryPrep.getClass().toString(), sampleId, e);
            System.exit(1);
        }
    }
}