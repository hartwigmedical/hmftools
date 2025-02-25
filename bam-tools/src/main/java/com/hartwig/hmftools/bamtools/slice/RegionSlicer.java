package com.hartwig.hmftools.bamtools.slice;

import static java.lang.String.format;

import static com.hartwig.hmftools.bamtools.common.CommonUtils.APP_NAME;
import static com.hartwig.hmftools.bamtools.common.CommonUtils.BT_LOGGER;
import static com.hartwig.hmftools.bamtools.slice.SliceConfig.UNMAPPED_READS_DISABLED;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource.loadRefGenome;
import static com.hartwig.hmftools.common.utils.PerformanceCounter.runTimeMinsStr;

import java.io.File;
import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.CompletableFuture;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.ThreadFactory;

import com.google.common.util.concurrent.ThreadFactoryBuilder;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;
import com.hartwig.hmftools.common.region.ChrBaseRegion;

import org.jetbrains.annotations.NotNull;

import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;


public class RegionSlicer
{
    private final SliceConfig mConfig;

    public RegionSlicer(final ConfigBuilder configBuilder)
    {
        mConfig = new SliceConfig(configBuilder);
    }

    public void run() throws ExecutionException, InterruptedException
    {
        if(!mConfig.isValid())
            System.exit(1);

        BT_LOGGER.info("starting BamSlicer");

        long startTimeMs = System.currentTimeMillis();

        SliceWriter sliceWriter = new SliceWriter(mConfig);
        ReadCache readCache = new ReadCache(sliceWriter);
        final ThreadLocal<SamReader> threadBamReader = createThreadLocalBamReader();
        final ExecutorService executorService = createExecutorService();

        final List<CompletableFuture<Void>> futures = new ArrayList<>();

        for(ChrBaseRegion region : mConfig.SliceRegions.Regions)
        {
            futures.add(CompletableFuture.runAsync(new RegionBamSlicer(region, mConfig, readCache, threadBamReader), executorService));
        }

        BT_LOGGER.info("splitting {} regions across {} threads", mConfig.SliceRegions.Regions.size(), mConfig.Threads);

        // wait for completion
        CompletableFuture.allOf(futures.toArray(CompletableFuture[]::new)).get();
        futures.clear();

        BT_LOGGER.info("initial slice complete");

        readCache.setProcessingRemoteRegions(true);

        // we need to perform remote position slicing twice. The reason is that we might still be missing supplementary reads of
        // mate reads that we get from the remote slicing
        for(int i = 0; i < 2; ++i)
        {
            List<ChrBaseRegion> remotePositions = readCache.getRemoteReadRegions();
            for(ChrBaseRegion region : remotePositions)
            {
                futures.add(CompletableFuture.runAsync(new RemoteReadSlicer(region, mConfig, readCache, threadBamReader), executorService));
            }
            BT_LOGGER.info("splitting {} remote regions across {} threads", remotePositions.size(), mConfig.Threads);

            // wait for completion
            CompletableFuture.allOf(futures.toArray(CompletableFuture[]::new)).get();
            futures.clear();
        }
        BT_LOGGER.info("remote slice complete");

        if(mConfig.MaxUnmappedReads != UNMAPPED_READS_DISABLED)
        {
            BT_LOGGER.info("slicing unmapped reads");

            sliceUnmappedReads(sliceWriter);

            BT_LOGGER.info("unmapped read slice complete");
        }

        sliceWriter.close();
        executorService.shutdown();

        if(mConfig.LogMissingReads)
        {
            readCache.logMissingReads();
        }

        BT_LOGGER.info("Regions slice complete, mins({})", runTimeMinsStr(startTimeMs));
    }

    private void sliceUnmappedReads(final SliceWriter sliceWriter)
    {
        SamReader samReader = !mConfig.RefGenomeFile.isEmpty() ?
                SamReaderFactory.makeDefault().referenceSequence(new File(mConfig.RefGenomeFile)).open(new File(mConfig.BamFile)) : null;

        long unmappedCount = 0;

        SAMRecordIterator iterator = samReader.queryUnmapped();

        while(iterator.hasNext())
        {
            sliceWriter.writeRead(iterator.next());
            ++unmappedCount;

            if(mConfig.MaxUnmappedReads > 0 && unmappedCount >= mConfig.MaxUnmappedReads)
                break;
        }
    }

    @NotNull
    private ExecutorService createExecutorService()
    {
        int numDigits = Integer.toString(mConfig.Threads - 1).length();
        final ThreadFactory namedThreadFactory = new ThreadFactoryBuilder().setNameFormat("thread-%0" + numDigits + "d").build();
        return Executors.newFixedThreadPool(mConfig.Threads, namedThreadFactory);
    }

    @NotNull
    private ThreadLocal<SamReader> createThreadLocalBamReader()
    {
        return ThreadLocal.withInitial(() ->
        {
            SamReaderFactory samReaderFactory = SamReaderFactory.makeDefault().validationStringency(ValidationStringency.LENIENT);
            if(mConfig.RefGenomeFile != null)
            {
                samReaderFactory = samReaderFactory.referenceSequence(new File(mConfig.RefGenomeFile));
            }
            return samReaderFactory.open(new File(mConfig.BamFile));
        });
    }

    public static void main(@NotNull final String[] args) throws ExecutionException, InterruptedException
    {
        ConfigBuilder configBuilder = new ConfigBuilder(APP_NAME);
        SliceConfig.addConfig(configBuilder);

        configBuilder.checkAndParseCommandLine(args);

        RegionSlicer regionSlicer = new RegionSlicer(configBuilder);
        regionSlicer.run();
    }
}
