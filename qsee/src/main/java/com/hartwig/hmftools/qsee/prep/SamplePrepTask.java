package com.hartwig.hmftools.qsee.prep;

import static com.hartwig.hmftools.qsee.common.QseeConstants.QC_LOGGER;

import java.nio.file.NoSuchFileException;
import java.util.ArrayList;
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
    private final boolean mAllowMissingInput;

    private final String mLogPrefix;

    @Nullable
    private final FeatureMatrix mSampleFeatureMatrix;

    public SamplePrepTask(CategoryPrep categoryPrep,
            List<String> sampleIds, int sampleIndex, SampleType sampleType,
            @Nullable FeatureMatrix sampleFeatureMatrix, boolean allowMissingInput)
    {
        mCategoryPrep = categoryPrep;
        mSampleIds = sampleIds;
        mSampleIndex = sampleIndex;
        mSampleType = sampleType;
        mSampleFeatureMatrix = sampleFeatureMatrix;
        mAllowMissingInput = allowMissingInput;

        mLogPrefix = logPrefix(sampleType, categoryPrep);
    }

    public static String logPrefix(SampleType sampleType, CategoryPrep categoryPrep)
    {
        return String.format("sampleType(%s) category(%s) -", sampleType, categoryPrep.name());
    }

    private void logProgress(int sampleIndex)
    {
        int sampleCount = mSampleIds.size();

        if(sampleCount == 1)
            return;

        int PROGRESS_INTERVAL = 100;

        boolean hasManySamples = sampleCount >= PROGRESS_INTERVAL*2;

        if(hasManySamples)
        {
            boolean isSampleAtInterval = (sampleIndex+1) % PROGRESS_INTERVAL == 0;
            boolean isLastSample = sampleIndex == sampleCount-1;

            if(isSampleAtInterval || isLastSample)
            {
                QC_LOGGER.debug("{} Progress: {}/{} - current sample: {}",
                        mLogPrefix, sampleIndex+1, sampleCount, mSampleIds.get(sampleIndex));
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
        }
        catch(NoSuchFileException e)
        {
            mFeatures = new ArrayList<>();
            QC_LOGGER.error("{} sample({}) missing input file(s): {}", mLogPrefix, sampleId, e.getMessage());

            if(!mAllowMissingInput)
                System.exit(1);
        }
        catch(Exception e)
        {
            QC_LOGGER.error("{} Failed to run prep for sample({})", mLogPrefix, sampleId, e);
            System.exit(1);
        }

        if(mSampleFeatureMatrix != null)
        {
            mSampleFeatureMatrix.addRow(sampleId, mFeatures);
            mFeatures = null;
        }
    }
}