package com.hartwig.hmftools.cup.prep;

import static com.hartwig.hmftools.cup.CuppaConfig.CUP_LOGGER;

import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.Callable;
import java.util.concurrent.ConcurrentHashMap;

import org.jetbrains.annotations.Nullable;

public class SampleOneCategoryTask implements Callable
{
    public final PrepConfig mConfig;
    public final CategoryPrep mCategoryPrep;
    public final int mSampleIndex;
    public final String mSampleName;

    @Nullable private List<DataItem> mDataItems;
    @Nullable ConcurrentHashMap<DataItem.Index, String[]> FeatureBySampleMatrix;

    public SampleOneCategoryTask(
            final int sampleIndex,
            final PrepConfig prepConfig,
            final CategoryPrep categoryPrep,
            @Nullable ConcurrentHashMap<DataItem.Index, String[]> featureBySampleMatrix
    ){
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

        if(mConfig.isMultiSample() & (totalSamples < 10 || sampleNum % 100 == 0))
            CUP_LOGGER.info("{}/{}: sample({})", sampleNum, totalSamples, mSampleName);

        mDataItems = mCategoryPrep.extractSampleData(mSampleName);
    }

    public List<DataItem> getDataItems()
    {
        return mDataItems;
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

    public void run()
    {
        processSample();

        if(mDataItems == null)
        {
            String errorMessage = String.format("Failed feature extraction category(%s) for sample(%s)", mCategoryPrep.categoryType(), mSampleName);

            if(mConfig.isSingleSample())
            {
                CUP_LOGGER.error(errorMessage);
                System.exit(1);
            }
            else
            {
                CUP_LOGGER.error(errorMessage + ". Output feature matrix will contain nulls for this sample");
            }
        }

        if(mDataItems != null & mConfig.isMultiSample())
        {
            addDataItemsToMatrix();
        }
    }

    @Override
    public Long call() throws Exception
    {
        run();
        return (long) 0;
    }
}
