package com.hartwig.hmftools.cup.prep;

import static com.hartwig.hmftools.cup.CuppaConfig.CUP_LOGGER;

import java.util.List;
import java.util.concurrent.Callable;
import java.util.concurrent.ConcurrentHashMap;

import org.jetbrains.annotations.Nullable;

public class SampleOneCategoryTask implements Callable
{
    public final PrepConfig mConfig;
    public final CategoryPrep mCategoryPrep;
    public final int mSampleIndex;

    private List<DataItem> mDataItems;
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

        if(mConfig.isMultiSample() & featureBySampleMatrix == null)
        {
            CUP_LOGGER.error("`featureBySampleMatrix` must not be null in multi sample mode");
            System.exit(1);
        }

        FeatureBySampleMatrix = featureBySampleMatrix;
    }

    public String getSampleId()
    {
        return mConfig.SampleIds.get(mSampleIndex);
    }

    public List<DataItem> processSample()
    {
        try
        {
            int sampleNum = mSampleIndex + 1;
            if(mConfig.isMultiSample() & (sampleNum % 100 == 0))
            {
                CUP_LOGGER.info("{}/{}: sample({})", sampleNum, mConfig.SampleIds.size(), getSampleId());
            }
            mDataItems = mCategoryPrep.extractSampleData(getSampleId());
        }
        catch(Exception e)
        {
            CUP_LOGGER.error("Feature extraction failed for category({})", mCategoryPrep.categoryType());
            System.exit(1);
        }

        return mDataItems;
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
            DataItem.Index featureIndex = dataItem.Index;

            FeatureBySampleMatrix.computeIfAbsent(featureIndex, k -> new String[nSamples]);
            FeatureBySampleMatrix.get(featureIndex)[mSampleIndex] = dataItem.Value;
        }
    }

    public void run()
    {
        processSample();

        if(mConfig.isMultiSample())
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
