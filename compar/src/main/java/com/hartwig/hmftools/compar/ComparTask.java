package com.hartwig.hmftools.compar;

import static com.hartwig.hmftools.compar.CommonUtils.buildComparers;
import static com.hartwig.hmftools.compar.Compar.writeSampleMismatches;
import static com.hartwig.hmftools.compar.ComparConfig.CMP_LOGGER;

import java.io.BufferedWriter;
import java.util.List;
import java.util.concurrent.Callable;

import com.google.common.collect.Lists;

public class ComparTask implements Callable
{
    private final int mTaskId;
    private final ComparConfig mConfig;
    private final List<String> mSampleIds;
    private final List<ItemComparer> mComparators;
    private final BufferedWriter mDiffWriter;

    public ComparTask(int taskId, final ComparConfig config, final BufferedWriter diffWriter)
    {
        mTaskId = taskId;
        mConfig = config;
        mDiffWriter = diffWriter;

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

            if(i > 0 && (i % 100) == 0)
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
        final List<Mismatch> mismatches = Lists.newArrayList();

        for(ItemComparer comparator : mComparators)
        {
            CMP_LOGGER.debug("sample({}) checking {}", sampleId, comparator.category());
            comparator.processSample(sampleId, mismatches);
        }

        CMP_LOGGER.debug("sample({}) writing {} mismatches", sampleId, mismatches.size());

        writeSampleMismatches(mDiffWriter, mConfig.multiSample() ? sampleId : null, mismatches);
    }

}
