package com.hartwig.hmftools.bamtools.synthetic;

import static java.lang.Math.max;
import static java.lang.Math.min;
import static java.lang.Math.round;

import static com.hartwig.hmftools.bamtools.common.CommonUtils.APP_NAME;
import static com.hartwig.hmftools.bamtools.common.CommonUtils.BT_LOGGER;
import static com.hartwig.hmftools.bamtools.synthetic.SyntheticConstants.MAX_CN_MULTIPLE;
import static com.hartwig.hmftools.common.perf.PerformanceCounter.runTimeMinsStr;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.concurrent.CompletableFuture;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.ThreadFactory;
import java.util.stream.Collectors;

import com.google.common.util.concurrent.ThreadFactoryBuilder;
import com.hartwig.hmftools.bamtools.slice.ReadCache;
import com.hartwig.hmftools.bamtools.slice.RegionBamSlicer;
import com.hartwig.hmftools.bamtools.slice.RemoteReadSlicer;
import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.region.ChrBaseRegion;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;

import org.jetbrains.annotations.NotNull;

import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;

public class SyntheticBamCreator
{
    private final SyntheticConfig mConfig;
    private final List<SamReader> mBamReaders;

    public SyntheticBamCreator(final ConfigBuilder configBuilder)
    {
        mConfig = new SyntheticConfig(configBuilder);
        mBamReaders = Collections.synchronizedList(new ArrayList<>());
    }

    public void run() throws ExecutionException, InterruptedException, IOException
    {
        BT_LOGGER.info("starting synthetic BAM creator");

        long startTimeMs = System.currentTimeMillis();

        BamWriter bamWriter = new BamWriter(mConfig);
        ReadCache readCache = new ReadCache(bamWriter, mConfig.Params.isDownsampling());

        RegionsBuilder regionsBuilder = new RegionsBuilder(mConfig);

        regionsBuilder.buildRegions();

        if(!mConfig.WriteBam)
        {
            BT_LOGGER.info("no BAM writing");
            return;
        }

        List<RegionData> regionDataList = regionsBuilder.regions();
        List<ChrBaseRegion> sliceRegions = regionDataList.stream().collect(Collectors.toList());

        ThreadLocal<SamReader> threadBamReader = createThreadLocalBamReader();
        ExecutorService executorService = createExecutorService();

        List<CompletableFuture<Void>> futures = new ArrayList<>();

        for(RegionData regionData : regionDataList)
        {
            RegionBamSlicer regionBamSlicer = new RegionBamSlicer(regionData, sliceRegions, mConfig.Params, readCache, threadBamReader);
            int regionTargetDepth = calculateTargetDepth(regionData);
            regionBamSlicer.setTargetDepth(regionTargetDepth);
            futures.add(CompletableFuture.runAsync(regionBamSlicer, executorService));
        }

        BT_LOGGER.info("splitting {} regions across {} threads", regionsBuilder.regions().size(), mConfig.Threads);

        // wait for completion
        CompletableFuture.allOf(futures.toArray(CompletableFuture[]::new)).get();
        futures.clear();

        BT_LOGGER.info("initial slice complete, read written({}) cached fragments({})",
                bamWriter.writeCount(), readCache.fragmentMap().size());

        readCache.setProcessingRemoteRegions(true);

        // perform the slice of remote positions twice if necessary to pick up remote mates which in turn have remote supplementaries
        for(int i = 0; i < 2; ++i)
        {
            List<ChrBaseRegion> remotePositions = readCache.collateRemoteReadRegions();

            for(ChrBaseRegion region : remotePositions)
            {
                if(!HumanChromosome.contains(region.Chromosome))
                    continue;

                futures.add(CompletableFuture.runAsync(
                        new RemoteReadSlicer(region, mConfig.Params, readCache, threadBamReader), executorService));
            }

            BT_LOGGER.info("splitting {} remote regions across {} threads", remotePositions.size(), mConfig.Threads);

            // wait for completion
            CompletableFuture.allOf(futures.toArray(CompletableFuture[]::new)).get();
            futures.clear();
        }

        BT_LOGGER.info("remote slice complete");

        executorService.shutdown();

        readCache.logMissingReads(100);

        readCache.flushCompleteFragments();

        bamWriter.close();
        closeBamReaders();

        BT_LOGGER.info("Synthetic BAM creation complete, mins({})", runTimeMinsStr(startTimeMs));
    }

    private int calculateTargetDepth(final RegionData regionData)
    {
        if(mConfig.Params.TargetDepth <= 0)
            return 0;

        int regionTargetDepth;

        if(mConfig.PurpleDir != null)
        {
            int expectedCopies = regionData.Chromosome.contains("Y") ? 1 : 2;
            double copyNumberAdjustment = min(max(0, regionData.copyNumber() / expectedCopies), MAX_CN_MULTIPLE);
            regionTargetDepth = (int)round(mConfig.Params.TargetDepth * copyNumberAdjustment);
        }
        else
        {
            regionTargetDepth = mConfig.Params.TargetDepth;
        }

        return regionTargetDepth;
    }

    private ExecutorService createExecutorService()
    {
        int numDigits = Integer.toString(mConfig.Threads - 1).length();
        final ThreadFactory namedThreadFactory = new ThreadFactoryBuilder().setNameFormat("thread-%0" + numDigits + "d").build();
        return Executors.newFixedThreadPool(mConfig.Threads, namedThreadFactory);
    }

    private ThreadLocal<SamReader> createThreadLocalBamReader()
    {
        return ThreadLocal.withInitial(() ->
        {
            SamReaderFactory samReaderFactory = SamReaderFactory.makeDefault().validationStringency(ValidationStringency.SILENT);
            if(mConfig.RefGenomeFile != null)
            {
                samReaderFactory = samReaderFactory.referenceSequence(new File(mConfig.RefGenomeFile));
            }
            SamReader bamReader = samReaderFactory.open(new File(mConfig.BamFile));
            mBamReaders.add(bamReader);
            return bamReader;
        });
    }

    private void closeBamReaders() throws IOException
    {
        for(SamReader bamReader : mBamReaders)
        {
            bamReader.close();
        }
    }

    public static void main(@NotNull final String[] args) throws ExecutionException, InterruptedException, IOException
    {
        ConfigBuilder configBuilder = new ConfigBuilder(APP_NAME);
        SyntheticConfig.addConfig(configBuilder);

        configBuilder.checkAndParseCommandLine(args);

        SyntheticBamCreator regionSlicer = new SyntheticBamCreator(configBuilder);
        regionSlicer.run();
    }
}
