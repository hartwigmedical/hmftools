package com.hartwig.hmftools.dnds.builder;

import static com.hartwig.hmftools.dnds.DndsCommon.DN_LOGGER;

import java.io.BufferedWriter;
import java.util.List;
import java.util.concurrent.Callable;
import java.util.concurrent.atomic.AtomicInteger;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.dnds.SampleMutationalLoad;
import com.hartwig.hmftools.dnds.SomaticVariant;

public class SampleTask implements Callable<Void>
{
    private final int mTaskId;
    private final List<String> mSampleIds;
    private final int mTotalSampleCount;
    private final int mThreads;
    private final AtomicInteger mProcessedCount;
    private final SampleDataLoader mSampleDataLoader;
    private final BufferedWriter mMutLoadWriter;
    private final BufferedWriter mVariantsWriter;

    public SampleTask(
            final int taskId, final int totalSampleCount, final int threads, final AtomicInteger processedCount,
            final SampleDataLoader sampleDataLoader, final BufferedWriter mutLoadWriter, final BufferedWriter variantsWriter)
    {
        mTaskId = taskId;
        mSampleIds = Lists.newArrayList();
        mTotalSampleCount = totalSampleCount;
        mThreads = threads;
        mProcessedCount = processedCount;
        mSampleDataLoader = sampleDataLoader;
        mMutLoadWriter = mutLoadWriter;
        mVariantsWriter = variantsWriter;
    }

    public void addSample(final String sampleId) { mSampleIds.add(sampleId); }
    public void addSamples(final List<String> sampleIds) { mSampleIds.addAll(sampleIds); }

    @Override
    public Void call()
    {
        for(String sampleId : mSampleIds)
        {
            processSample(sampleId);

            int processed = mProcessedCount.incrementAndGet();

            if(mTotalSampleCount < 10 || (processed % 10) == 0)
            {
                DN_LOGGER.info("processed {} of {} samples", processed, mTotalSampleCount);
            }
        }

        if(mThreads > 1)
        {
            DN_LOGGER.debug("{}: task complete for {} samples", mTaskId, mSampleIds.size());
        }

        return null;
    }

    private void processSample(final String sampleId)
    {
        SampleData sampleData = mSampleDataLoader.loadSampleData(sampleId);

        if(sampleData == null)
            return;

        DN_LOGGER.debug("sample({}) loaded {} variants", sampleId, sampleData.Variants.size());

        SampleMutationalLoad.writeSampleMutationalLoad(mMutLoadWriter, sampleId, sampleData.MutationalLoad);
        SomaticVariant.writeVariants(mVariantsWriter, sampleId, sampleData.Variants);
    }
}
