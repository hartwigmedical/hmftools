package com.hartwig.hmftools.qsee.prep;

import static com.hartwig.hmftools.qsee.common.QseeConstants.QC_LOGGER;

import java.nio.file.NoSuchFileException;
import java.util.ArrayList;
import java.util.List;

import org.jetbrains.annotations.Nullable;

import com.hartwig.hmftools.qsee.cohort.FeatureMatrix;
import com.hartwig.hmftools.qsee.common.SampleType;
import com.hartwig.hmftools.qsee.feature.Feature;

public class CategoryPrepTask implements Runnable
{
    private final CategoryPrep mCategoryPrep;

    private final String mSampleId;
    private final int mSampleIndex;
    private final int mTotalSampleCount;
    private final SampleType mSampleType;
    private final boolean mAllowMissingInput;

    private final String mLogPrefix;

    private List<Feature> mOutput;

    @Nullable
    private final FeatureMatrix mSampleFeatureMatrix;

    public CategoryPrepTask(CategoryPrep categoryPrep,
            String sampleId, int sampleIndex, int totalSampleCount, SampleType sampleType,
            @Nullable FeatureMatrix sampleFeatureMatrix, boolean allowMissingInput)
    {
        mCategoryPrep = categoryPrep;
        mSampleId = sampleId;
        mSampleIndex = sampleIndex;
        mTotalSampleCount = totalSampleCount;
        mSampleType = sampleType;
        mSampleFeatureMatrix = sampleFeatureMatrix;
        mAllowMissingInput = allowMissingInput;

        mLogPrefix = logPrefix(sampleType, categoryPrep);
    }

    public static String logPrefix(SampleType sampleType, CategoryPrep categoryPrep)
    {
        return String.format("sampleType(%s) category(%s) -", sampleType, categoryPrep.name());
    }

    private void logProgress()
    {
        if(mTotalSampleCount == 1)
            return;

        int PROGRESS_INTERVAL = 100;

        boolean hasManySamples = mTotalSampleCount >= PROGRESS_INTERVAL * 2;

        if(hasManySamples)
        {
            boolean isSampleAtInterval = (mSampleIndex + 1) % PROGRESS_INTERVAL == 0;
            boolean isLastSample = mSampleIndex == mTotalSampleCount - 1;

            if(isSampleAtInterval || isLastSample)
            {
                QC_LOGGER.debug("{} Progress: {}/{} - current sample: {}",
                        mLogPrefix, mSampleIndex + 1, mTotalSampleCount, mSampleId);
            }
        }
    }

    @Override
    public void run()
    {
        try
        {
            logProgress();
            mOutput = mCategoryPrep.extractSampleData(mSampleId, mSampleType);
        }
        catch(NoSuchFileException e)
        {
            QC_LOGGER.error("{} sample({}) missing input file(s): {}", mLogPrefix, mSampleId, e.getMessage());

            if(!mAllowMissingInput)
                System.exit(1);

            mOutput = new ArrayList<>();
        }
        catch(Exception e)
        {
            QC_LOGGER.error("{} Failed to run prep for sample({})", mLogPrefix, mSampleId, e);
            System.exit(1);
        }

        if(mSampleFeatureMatrix != null)
        {
            mSampleFeatureMatrix.addRow(mSampleId, mOutput);
            mOutput = null;
        }
    }

    public List<Feature> getOutput(){ return mOutput; }
}