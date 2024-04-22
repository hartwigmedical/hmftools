package com.hartwig.hmftools.bamtools.btofmc;

import static java.lang.Math.max;

import static com.hartwig.hmftools.common.region.PartitionUtils.partitionChromosome;

import java.io.File;
import java.nio.charset.StandardCharsets;
import java.util.List;
import java.util.Queue;
import java.util.concurrent.ConcurrentLinkedQueue;
import java.util.concurrent.atomic.AtomicInteger;
import java.util.function.BiConsumer;
import java.util.function.Consumer;
import java.util.function.Supplier;

import com.google.common.io.CharSink;
import com.google.common.io.Files;
import com.google.inject.AbstractModule;
import com.google.inject.Provides;
import com.google.inject.Singleton;
import com.google.inject.TypeLiteral;
import com.google.inject.util.Providers;
import com.hartwig.hmftools.bamtools.btofmc.BindingAnnotations.Fastq1;
import com.hartwig.hmftools.bamtools.btofmc.BindingAnnotations.Fastq2;
import com.hartwig.hmftools.bamtools.btofmc.BindingAnnotations.TotalRegionCount;
import com.hartwig.hmftools.bamtools.btofmc.readcache.ConcurrentPartitionedReadCache;
import com.hartwig.hmftools.bamtools.btofmc.readcache.LockableReadCache;
import com.hartwig.hmftools.bamtools.btofmc.readcache.ReadCache;
import com.hartwig.hmftools.bamtools.btofmc.readcache.ReadCacheInterface;
import com.hartwig.hmftools.bamtools.btofmc.writer.ConcurrentPairedFastqWriter;
import com.hartwig.hmftools.bamtools.btofmc.writer.NoSyncPairedFastqWriter;
import com.hartwig.hmftools.bamtools.btofmc.writer.PairedFastqWriterInterface;
import com.hartwig.hmftools.common.region.ChrBaseRegion;

import org.apache.commons.compress.utils.Lists;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;

// TODO NEXT: TEST
// TODO: No logic in modules.
// TODO: Guice trace object creation.
// TODO: Listener for closing a flushing singletons, and check that closable non singleton instances are not provided.
public class DefaultModule extends AbstractModule
{
    private final BamToFastqConfig mConfig;
    private final Queue<RegionTask> mRegionTasks;
    private final int mTotalRegionTaskCount;

    public DefaultModule(final BamToFastqConfig config)
    {
        mConfig = config;
        mRegionTasks = getRegionTasks(config);
        mTotalRegionTaskCount = mRegionTasks.size();
    }

    private static Queue<RegionTask> getRegionTasks(final BamToFastqConfig config)
    {
        Queue<RegionTask> regionTasks = new ConcurrentLinkedQueue<>();
        regionTasks.add(RegionTask.UNMAPPED_READS);
        for(SAMSequenceRecord refSequence : config.RefGenome.refGenomeFile().getSequenceDictionary().getSequences())
        {
            for(ChrBaseRegion partition : partitionChromosome(refSequence, config.PartitionSize))
            {
                regionTasks.add(RegionTask.createChrRegion(partition));
            }
        }

        return regionTasks;
    }

    @Provides
    @Singleton
    @TotalRegionCount
    public int provideTotalRegionTaskCount()
    {
        return mTotalRegionTaskCount;
    }

    // TODO: clean this up.
    @Provides
    @Singleton
    public Supplier<LockableReadCache> provideLockableReadCacheSupplier(final BiConsumer<SAMRecord, SAMRecord> readPairConsumer)
    {
        return () -> new LockableReadCache(mConfig, () -> new ReadCache(mConfig, readPairConsumer));
    }

    @Provides
    @Singleton
    public List<RegionTaskConsumer> provideRegionTaskConsumers(final RegionTaskConsumerFactory regionTaskConsumerFactory)
    {
        List<RegionTaskConsumer> regionTaskConsumers = Lists.newArrayList();
        for(int i = 0; i < max(1, mConfig.Threads); ++i)
        {
            regionTaskConsumers.add(regionTaskConsumerFactory.create());
        }

        return regionTaskConsumers;
    }

    @Override
    public void configure()
    {
        bind(BamToFastqConfig.class).toInstance(mConfig);
        bind(new TypeLiteral<Queue<RegionTask>>() {}).toInstance(mRegionTasks);
        bind(new TypeLiteral<Supplier<SamReader>>() {}).toInstance(
                () -> SamReaderFactory.makeDefault().referenceSequence(new File(mConfig.RefGenomeFile)).open(new File(mConfig.BamFile)));
        bind(AtomicInteger.class).toInstance(new AtomicInteger());

        bind(RegionTaskConsumerFactory.class).in(Singleton.class);

        if(mConfig.Threads <= 1)
        {
            if(mConfig.NoWrite)
            {
                bind(CharSink.class).annotatedWith(Fastq1.class).toProvider(Providers.of(null));
                bind(CharSink.class).annotatedWith(Fastq2.class).toProvider(Providers.of(null));
            }
            else
            {
                CharSink fastq1 = Files.asCharSink(new File(mConfig.Fastq1OutputFile), StandardCharsets.UTF_8);
                CharSink fastq2 = Files.asCharSink(new File(mConfig.Fastq2OutputFile), StandardCharsets.UTF_8);
                bind(CharSink.class).annotatedWith(Fastq1.class).toInstance(fastq1);
                bind(CharSink.class).annotatedWith(Fastq2.class).toInstance(fastq2);
            }

            bind(ReadCache.class).in(Singleton.class);
            bind(ReadCacheInterface.class).to(ReadCache.class).in(Singleton.class);

            bind(NoSyncPairedFastqWriter.class).in(Singleton.class);
            bind(PairedFastqWriterInterface.class).to(NoSyncPairedFastqWriter.class).in(Singleton.class);
        }
        else
        {
            bind(ConcurrentPartitionedReadCache.class).in(Singleton.class);
            bind(ReadCacheInterface.class).to(ConcurrentPartitionedReadCache.class).in(Singleton.class);

            bind(ConcurrentPairedFastqWriter.class).in(Singleton.class);
            bind(PairedFastqWriterInterface.class).to(ConcurrentPairedFastqWriter.class).in(Singleton.class);
        }

        bind(new TypeLiteral<Consumer<SAMRecord>>() {}).to(ReadCacheInterface.class).in(Singleton.class);

        bind(new TypeLiteral<BiConsumer<SAMRecord, SAMRecord>>() {}).to(PairedFastqWriterInterface.class).in(Singleton.class);
    }
}
