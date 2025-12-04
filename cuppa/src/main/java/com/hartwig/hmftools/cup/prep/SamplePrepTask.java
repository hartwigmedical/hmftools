package com.hartwig.hmftools.cup.prep;

import static com.hartwig.hmftools.cup.common.CupConstants.CUP_LOGGER;

import java.util.List;
import java.util.concurrent.Callable;
import java.util.concurrent.ConcurrentHashMap;

import org.jetbrains.annotations.Nullable;

public class SamplePrepTask implements Callable<Void>
{
    private final PrepConfig mConfig;
    private final int mSampleIndex;
    private final String mSampleId;

    private CategoryPrep mCategoryPrep;

    @Nullable private List<DataItem> mDataItems;
    @Nullable private ConcurrentHashMap<DataItem.Index, String[]> FeatureBySampleMatrix;

    private static final int PROGRESS_INTERVAL = 100;

    public SamplePrepTask(
            final int sampleIndex,
            final PrepConfig prepConfig,
            final CategoryPrep categoryPrep,
            @Nullable final ConcurrentHashMap<DataItem.Index, String[]> featureBySampleMatrix)
    {
        mConfig = prepConfig;
        mCategoryPrep = categoryPrep;
        mSampleIndex = sampleIndex;

        mSampleId = categoryPrep.categoryType().isRna()
                ? mConfig.RnaSampleIds.get(mSampleIndex)
                : mConfig.SampleIds.get(mSampleIndex);

        if(mConfig.isMultiSample() & featureBySampleMatrix == null)
        {
            CUP_LOGGER.error("`featureBySampleMatrix` must not be null in multi sample mode");
            System.exit(1);
        }

        FeatureBySampleMatrix = featureBySampleMatrix;
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

    public List<DataItem> dataItems() { return mDataItems; }

    public void clearDataItems() { mDataItems = null; }

    public void clearCategoryPrep() { mCategoryPrep = null; }

    public void run()
    {
        if(mConfig.isMultiSample() & (mSampleIndex < PROGRESS_INTERVAL || mSampleIndex % PROGRESS_INTERVAL == 0))
        {
            int sampleNum = mSampleIndex + 1;
            int totalSamples = mConfig.SampleIds.size();
            CUP_LOGGER.debug("{}/{}: sample({})", sampleNum, totalSamples, mSampleId);
        }

        boolean skipRnaCategory = mCategoryPrep.categoryType().isRna() && mSampleId.equals(SampleIdsLoader.NO_RNA_SAMPLE_ID);
        if(!skipRnaCategory)
            mDataItems = mCategoryPrep.extractSampleData(mSampleId);

        if(mConfig.isMultiSample())
        {
            if(mDataItems == null)
            {
                // The matrix header will always use `SampleId`, and not `RnaSampleId`
                String sampleIdForLog =  mConfig.SampleIds.get(mSampleIndex);

                CUP_LOGGER.warn("multi-sample feature matrix will contain nulls for sample({}) category({})",
                        sampleIdForLog, mCategoryPrep.categoryType());
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
    public Void call()
    {
        run();
        return null;
    }
}
