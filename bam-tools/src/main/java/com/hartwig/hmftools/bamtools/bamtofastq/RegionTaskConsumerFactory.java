package com.hartwig.hmftools.bamtools.bamtofastq;

import java.util.concurrent.atomic.AtomicInteger;
import java.util.function.Consumer;
import java.util.function.Supplier;

import com.google.inject.Inject;
import com.hartwig.hmftools.bamtools.bamtofastq.collection.AtomicPStack;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReader;

public class RegionTaskConsumerFactory
{
    private final BamToFastqConfig mConfig;
    private final Supplier<SamReader> mSamReaderSupplier;
    private final Consumer<SAMRecord> mReadCache;
    private final AtomicPStack<RegionTask> mRegionTasks;
    private final AtomicInteger mTotalRegionCount;

    @Inject
    public RegionTaskConsumerFactory(final BamToFastqConfig config, final Supplier<SamReader> samReaderSupplier,
            final Consumer<SAMRecord> readCache, final AtomicPStack<RegionTask> regionTasks,
            final AtomicInteger totalRegionCount)
    {
        mConfig = config;
        mSamReaderSupplier = samReaderSupplier;
        mReadCache = readCache;
        mRegionTasks = regionTasks;
        mTotalRegionCount = totalRegionCount;
    }

    public RegionTaskConsumer create()
    {
        RegionTaskConsumer consumer = new RegionTaskConsumer(mConfig, mSamReaderSupplier, mReadCache, mRegionTasks, mTotalRegionCount);
        consumer.setUncaughtExceptionHandler(BamToFastq::uncaughtExceptionHandler);
        return consumer;
    }
}
