package com.hartwig.hmftools.qsee.prep;

import static com.hartwig.hmftools.qsee.common.QseeConstants.QC_LOGGER;

import java.util.List;

import org.jetbrains.annotations.Nullable;

import com.hartwig.hmftools.qsee.cohort.FeatureMatrix;
import com.hartwig.hmftools.qsee.common.SampleType;
import com.hartwig.hmftools.qsee.feature.Feature;

public class SamplePrepTask implements Runnable
{
    private final CategoryPrep mCategoryPrep;

    private final List<String> mSampleIds;
    private final int mSampleIndex;

    private final SampleType mSampleType;

    private List<Feature> mFeatures;

    @Nullable
    private final FeatureMatrix mSampleFeatureMatrix;

    public SamplePrepTask(CategoryPrep categoryPrep,
            List<String> sampleIds, int sampleIndex, SampleType sampleType,
            @Nullable FeatureMatrix sampleFeatureMatrix)
    {
        mCategoryPrep = categoryPrep;

        mSampleIds = sampleIds;
        mSampleIndex = sampleIndex;
        mSampleType = sampleType;

        mSampleFeatureMatrix = sampleFeatureMatrix;
    }

    private void logProgress(int sampleIndex)
    {
        int sampleCount = mSampleIds.size();

        if(sampleCount == 1)
            return;

        int PROGRESS_INTERVAL = 100;

        boolean hasManySamples = sampleCount >= PROGRESS_INTERVAL;

        if(hasManySamples)
        {
            boolean isSampleAtInterval = (sampleIndex+1) % PROGRESS_INTERVAL == 0;
            boolean isLastSample = sampleIndex == sampleCount-1;

            if(isLastSample || isSampleAtInterval)
            {
                QC_LOGGER.debug("Progress: {}/{} - current sample: {}", sampleIndex+1, sampleCount, mSampleIds.get(sampleIndex));
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

            mFeatures = mCategoryPrep.extractSampleData(sampleId, mSampleType);

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