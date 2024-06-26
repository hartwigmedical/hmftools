package com.hartwig.hmftools.cup.prep;

import static com.hartwig.hmftools.cup.common.CupConstants.CUP_LOGGER;

import java.util.List;
import java.util.concurrent.Callable;
import java.util.concurrent.ConcurrentHashMap;

import org.jetbrains.annotations.Nullable;

public class SampleOneCategoryTask implements Callable
{
    public final PrepConfig mConfig;
    public final int mSampleIndex;
    public final String mSampleName;

    @Nullable public CategoryPrep mCategoryPrep;

    @Nullable public List<DataItem> mDataItems;
    @Nullable public ConcurrentHashMap<DataItem.Index, String[]> FeatureBySampleMatrix;

    public SampleOneCategoryTask(
            final int sampleIndex,
            final PrepConfig prepConfig,
            CategoryPrep categoryPrep,
            @Nullable ConcurrentHashMap<DataItem.Index, String[]> featureBySampleMatrix)
    {
        mConfig = prepConfig;
        mCategoryPrep = categoryPrep;
        mSampleIndex = sampleIndex;
        mSampleName = mConfig.SampleIds.get(mSampleIndex);

        if(mConfig.isMultiSample() & featureBySampleMatrix == null)
        {
            CUP_LOGGER.error("`featureBySampleMatrix` must not be null in multi sample mode");
            System.exit(1);
        }

        FeatureBySampleMatrix = featureBySampleMatrix;
    }

    public void processSample()
    {
        int sampleNum = mSampleIndex + 1;
        int totalSamples = mConfig.SampleIds.size();

        if(mConfig.isMultiSample() & sampleNum % mConfig.ProgressInterval == 0)
        {
            CUP_LOGGER.info("{}/{}: sample({})", sampleNum, totalSamples, mSampleName);
        }

        mDataItems = mCategoryPrep.extractSampleData(mSampleName);
    }

    public synchronized void addDataItemsToMatrix()
    {
        int nSamples = mConfig.SampleIds.size();

        for(DataItem dataItem : mDataItems)
        {
            FeatureBySampleMatrix.computeIfAbsent(dataItem.Index, k -> new String[nSamples]);
            FeatureBySampleMatrix.get(dataItem.Index)[mSampleIndex] = dataItem.Value;
        }
    }

    public void clearDataItems()
    {
        mDataItems = null;
    }

    public void clearCategoryPrep()
    {
        mCategoryPrep = null;
    }

    public void run()
    {
        processSample();

        if(mConfig.isMultiSample())
        {
            if(mDataItems == null)
            {
                CUP_LOGGER.error("multi-sample feature matrix will contain nulls for sample({}) category({})", mSampleName, mCategoryPrep.categoryType());
            }
            else
            {
                addDataItemsToMatrix();
                clearDataItems();
            }
        }

        clearCategoryPrep();
    }

    @Override
    public Long call()
    {
        run();
        return (long) 0;
    }
}
