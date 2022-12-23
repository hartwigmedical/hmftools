package com.hartwig.hmftools.compar;

import static java.lang.String.format;

import static com.hartwig.hmftools.compar.CommonUtils.buildComparers;
import static com.hartwig.hmftools.compar.ComparConfig.CMP_LOGGER;

import java.util.List;
import java.util.concurrent.Callable;

import com.google.common.collect.Lists;

public class ComparTask implements Callable
{
    private final int mTaskId;
    private final ComparConfig mConfig;
    private final List<String> mSampleIds;
    private final List<ItemComparer> mComparators;

    private final MismatchWriter mWriter;

    public ComparTask(int taskId, final ComparConfig config, final MismatchWriter writer)
    {
        mTaskId = taskId;
        mConfig = config;
        mWriter = writer;

        mSampleIds = Lists.newArrayList();
        mComparators = buildComparers(config);
    }

    public List<String> getSampleIds() { return mSampleIds; }

    @Override
    public Long call()
    {
        for(int i = 0; i < mSampleIds.size(); ++i)
        {
            String sampleId = mSampleIds.get(i);

            processSample(sampleId);

            if(i > 0 && (i % 10) == 0)
            {
                // clean-up memory esp from somatic variant DB retrieval
                System.gc();
            }

            if(i > 0 && (i % 50) == 0)
            {
                CMP_LOGGER.info("{}: processed {} samples", mTaskId, i);
            }
        }

        if(mConfig.Threads > 1)
        {
            CMP_LOGGER.info("{}: tasks complete for {} samples", mTaskId, mSampleIds.size());
        }

        return (long)0;
    }

    private void processSample(final String sampleId)
    {
        int totalMismatches = 0;
        int failedTypes = 0;
        for(ItemComparer comparer : mComparators)
        {
            List<Mismatch> mismatches = Lists.newArrayList();

            try
            {
                boolean status = comparer.processSample(sampleId, mismatches);

                if(!status)
                    ++failedTypes;
            }
            catch(Exception e)
            {
                CMP_LOGGER.error("sample({}) failed processing: {}", sampleId, e.toString());
                e.printStackTrace();
                ++failedTypes;
            }

            mWriter.writeSampleMismatches(sampleId, comparer, mismatches);
            totalMismatches += mismatches.size();
        }

        CMP_LOGGER.debug("sample({}) wrote {} mismatches {}",
                sampleId, totalMismatches, failedTypes > 0 ? format("failed types %d", failedTypes) : "");
    }
}
