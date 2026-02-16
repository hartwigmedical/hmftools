package com.hartwig.hmftools.qsee.prep;

import static com.hartwig.hmftools.qsee.common.QseeConstants.QC_LOGGER;

import java.io.IOException;
import java.nio.file.NoSuchFileException;
import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.atomic.AtomicInteger;

import org.apache.logging.log4j.Level;
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

    @Nullable
    private final AtomicInteger mSamplesMissingInputCount;

    private List<Feature> mOutput;

    @Nullable
    private final FeatureMatrix mSampleFeatureMatrix;

    public CategoryPrepTask(CategoryPrep categoryPrep,
            String sampleId, int sampleIndex, int totalSampleCount, SampleType sampleType,
            @Nullable FeatureMatrix sampleFeatureMatrix, boolean allowMissingInput, @Nullable AtomicInteger samplesMissingInputCount)
    {
        mCategoryPrep = categoryPrep;
        mSampleId = sampleId;
        mSampleIndex = sampleIndex;
        mTotalSampleCount = totalSampleCount;
        mSampleType = sampleType;
        mSampleFeatureMatrix = sampleFeatureMatrix;
        mAllowMissingInput = allowMissingInput;
        mSamplesMissingInputCount = samplesMissingInputCount;
    }

    public CategoryPrepTask(CategoryPrep categoryPrep, String sampleId, SampleType sampleType, boolean allowMissingInput)
    {
        mCategoryPrep = categoryPrep;
        mSampleId = sampleId;
        mSampleIndex = 0;
        mTotalSampleCount = 0;
        mSampleType = sampleType;
        mSampleFeatureMatrix = null;
        mAllowMissingInput = allowMissingInput;
        mSamplesMissingInputCount = null;
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
                QC_LOGGER.debug("category({}) - Progress: {}/{} - current sample: {}",
                        mCategoryPrep.name(), mSampleIndex + 1, mTotalSampleCount, mSampleId);
            }
        }
    }

    public static void missingInputFilesError(
            boolean allowMissingInput, CategoryPrep categoryPrep, SampleType sampleType, String sampleId, String missingFilePath)
    {
        QC_LOGGER.log(
                allowMissingInput ? Level.WARN : Level.ERROR,
                "sampleType({}) category({}) - sample({}) missing input file: {}",
                sampleType, categoryPrep.name(), sampleId, missingFilePath
        );

        if(!allowMissingInput)
            System.exit(1);
    }

    @Override
    public void run()
    {
        try
        {
            logProgress();
            mOutput = mCategoryPrep.extractSampleData(mSampleId, mSampleType);
        }
        catch(IOException e)
        {
            missingInputFilesError(mAllowMissingInput, mCategoryPrep, mSampleType, mSampleId, e.getMessage());

            if(mSamplesMissingInputCount != null)
                mSamplesMissingInputCount.incrementAndGet();

            mOutput = new ArrayList<>();
        }
        catch(Exception e)
        {
            QC_LOGGER.error("sampleType({}) category({}) - Failed to run prep for sample({})", mCategoryPrep.name(), mSampleId, e);
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