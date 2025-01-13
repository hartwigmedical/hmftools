package com.hartwig.hmftools.redux;

import static java.lang.Math.max;

import static com.hartwig.hmftools.common.utils.PerformanceCounter.runTimeMinsStr;
import static com.hartwig.hmftools.common.utils.TaskExecutor.runThreadTasks;
import static com.hartwig.hmftools.redux.PartitionThread.splitRegionsIntoPartitions;
import static com.hartwig.hmftools.redux.ReduxConfig.APP_NAME;
import static com.hartwig.hmftools.redux.ReduxConfig.RD_LOGGER;
import static com.hartwig.hmftools.redux.ReduxConfig.registerConfig;
import static com.hartwig.hmftools.redux.common.Constants.DEFAULT_READ_LENGTH;
import static com.hartwig.hmftools.redux.unmap.RegionUnmapper.createThreadTasks;
import static com.hartwig.hmftools.redux.write.PartitionInfo.partitionInfoStr;

import java.util.Collections;
import java.util.List;
import java.util.Queue;
import java.util.concurrent.ConcurrentLinkedQueue;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.bamops.BamSampler;
import com.hartwig.hmftools.common.basequal.jitter.ConsensusMarker;
import com.hartwig.hmftools.common.basequal.jitter.JitterAnalyser;
import com.hartwig.hmftools.common.basequal.jitter.JitterAnalyserConfig;
import com.hartwig.hmftools.common.region.ChrBaseRegion;
import com.hartwig.hmftools.common.utils.PerformanceCounter;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;
import com.hartwig.hmftools.redux.common.Statistics;
import com.hartwig.hmftools.redux.unmap.RegionUnmapper;
import com.hartwig.hmftools.redux.unmap.TaskQueue;
import com.hartwig.hmftools.redux.unmap.UnmapStats;
import com.hartwig.hmftools.redux.write.FileWriterCache;
import com.hartwig.hmftools.redux.write.FinalBamWriter;
import com.hartwig.hmftools.redux.write.PartitionInfo;

public class ReduxApplication
{
    private final ReduxConfig mConfig;

    public ReduxApplication(final ConfigBuilder configBuilder)
    {
        mConfig = new ReduxConfig(configBuilder);
    }

    public void run()
    {
        if(!mConfig.isValid())
            System.exit(1);

        RD_LOGGER.info("sample({}) starting duplicate marking", mConfig.SampleId);

        long startTimeMs = System.currentTimeMillis();

        setReadLength();

        JitterAnalyser jitterAnalyser = null;

        if(mConfig.JitterMsiFile != null)
        {
            JitterAnalyserConfig jitterConfig = new JitterAnalyserConfig(
                    mConfig.SampleId, mConfig.RefGenVersion, mConfig.RefGenomeFile, mConfig.JitterMsiFile, mConfig.OutputDir,
                    mConfig.JitterMsiMaxSitePercContribution, false);

            ConsensusMarker consensusMarker = ConsensusMarker.fromSequencingType(mConfig.Sequencing);
            jitterAnalyser = new JitterAnalyser(jitterConfig, RD_LOGGER, consensusMarker);
        }

        FileWriterCache fileWriterCache = new FileWriterCache(mConfig, jitterAnalyser);
        UnmapStats unmapStats = mConfig.UnmapRegions.stats();

        if(mConfig.UnmapRegions.enabled())
        {
            List<Thread> unmappingThreadTasks = Lists.newArrayList();
            List<RegionUnmapper> readUnmappers = createThreadTasks(mConfig, fileWriterCache, unmappingThreadTasks);

            if(!readUnmappers.isEmpty())
            {
                if(!runThreadTasks(unmappingThreadTasks))
                    System.exit(1);

                RD_LOGGER.debug("initial unmapping complete");

                long readsProcessed = readUnmappers.stream().mapToLong(x -> x.processedReads()).sum();
                RD_LOGGER.info("readsProcessed({}) unmapped stats: {}", readsProcessed, mConfig.UnmapRegions.stats());

                // reset unmapped stats for a final comparison
                mConfig.UnmapRegions.setStats(new UnmapStats());
            }

            if(!fileWriterCache.prepareSortedUnmappingBam())
                System.exit(1);
        }

        // partition the genome into sequential regions to be processed by each thread
        List<PartitionThread> partitionThreads = createPartitionThreads(fileWriterCache);
        List<Thread> allThreads = Lists.newArrayList(partitionThreads);

        FinalBamWriter finalBamWriter = null;

        if(mConfig.WriteBam && mConfig.ParallelConcatenation)
        {
            finalBamWriter = new FinalBamWriter(mConfig, fileWriterCache);
            allThreads.add(finalBamWriter);
        }

        if(!runThreadTasks(allThreads))
            System.exit(1);

        RD_LOGGER.info("all partition tasks complete");

        if(jitterAnalyser != null)
        {
            try
            {
                RD_LOGGER.info("analysing microsatellite jitter");
                jitterAnalyser.writeAnalysisOutput();
            }
            catch(Exception e)
            {
                RD_LOGGER.error("failed to write output of jitter analysis: {}", e.toString());
                System.exit(1);
            }
        }

        if(mConfig.JitterMsiOnly)
        {
            RD_LOGGER.info("Redux jitter complete, mins({})", runTimeMinsStr(startTimeMs));
            return;
        }

        List<PartitionReader> partitionReaders = partitionThreads.stream().map(x -> x.partitionReader()).collect(Collectors.toList());

        Statistics combinedStats = new Statistics();
        partitionReaders.forEach(x -> combinedStats.merge(x.statistics()));

        List<PerformanceCounter> combinedPerfCounters = mergePerfCounters(partitionReaders);

        // free up any processing state
        partitionReaders.clear();

        // usually avoid manual calls to this but since the external BAM tools make independent calls to access memory and
        // the core routines are complete, it is helpful to do so now
        System.gc();

        fileWriterCache.finaliseBams();

        long sortedBamUnsortedWriteCount = fileWriterCache.sortedBamUnsortedWriteCount();

        if(sortedBamUnsortedWriteCount > 0)
        {
            RD_LOGGER.warn("unsorted BAM write count({}) via sorted BAM writers", sortedBamUnsortedWriteCount);
        }

        combinedStats.logStats();

        if(mConfig.UnmapRegions.enabled())
        {
            if(mConfig.RunChecks)
                mConfig.readChecker().logUnmatchedUnmappedReads();

            // check that the unmapping counts match the re-tested unmapped reads from the the partition readers
            UnmapStats reunmapStats = mConfig.UnmapRegions.stats();

            RD_LOGGER.debug("re-unmapped stats: {}", reunmapStats.toString());

            if(reunmapStats.ReadCount.get() != unmapStats.ReadCount.get()
            || reunmapStats.FullyUnmappedCount.get() != unmapStats.FullyUnmappedCount.get())
            {
                RD_LOGGER.warn("re-unmapped stats differ: {}", reunmapStats.toString());
            }
        }

        if(mConfig.WriteStats)
        {
            combinedStats.writeDuplicateStats(mConfig);

            if(mConfig.UMIs.Enabled)
            {
                combinedStats.UmiStats.writePositionFragmentsData(mConfig);

                if(mConfig.UMIs.BaseStats)
                {
                    combinedStats.UmiStats.writeUmiBaseDiffStats(mConfig);
                    combinedStats.UmiStats.writeUmiBaseFrequencyStats(mConfig);
                }
            }
        }

        if(finalBamWriter != null)
            finalBamWriter.logTimes();

        logPerformanceStats(combinedPerfCounters);

        RD_LOGGER.info("Redux complete, mins({})", runTimeMinsStr(startTimeMs));
    }

    private List<PartitionThread> createPartitionThreads(final FileWriterCache fileWriterCache)
    {
        int partitionThreadCount = mConfig.WriteBam & mConfig.ParallelConcatenation ? max(mConfig.Threads - 1, 1) : mConfig.Threads;
        int partitionCount = mConfig.PartitionThreadRatio * partitionThreadCount;

        RD_LOGGER.debug("splitting {} partition regions across {} threads", partitionCount, partitionThreadCount);

        List<PartitionThread> partitionThreads = Lists.newArrayListWithCapacity(partitionThreadCount);

        List<List<ChrBaseRegion>> partitionRegions = splitRegionsIntoPartitions(
                mConfig.SpecificChrRegions, partitionCount, mConfig.RefGenVersion, mConfig.RefGenome);

        if(partitionRegions.isEmpty())
            return Collections.emptyList();

        for(List<ChrBaseRegion> regions : partitionRegions)
        {
            if(regions.isEmpty())
                break;

            long regionsLength = regions.stream().mapToLong(x -> x.baseLength()).sum();

            RD_LOGGER.debug("adding partition regions({}) totalLength({}): {}", regions.size(), regionsLength, partitionInfoStr(regions));

            fileWriterCache.addPartition(regions);
        }

        Queue<PartitionInfo> partitionQueue = new ConcurrentLinkedQueue<>();
        partitionQueue.addAll(fileWriterCache.partitions());

        TaskQueue taskQueue = new TaskQueue(partitionQueue, "partitions", 0); // log on completed partitions

        List<String> inputBamFiles = Lists.newArrayList(mConfig.BamFiles);

        if(mConfig.UnmapRegions.enabled() && mConfig.WriteBam)
            inputBamFiles.add(fileWriterCache.unmappedSortedBamFilename());

        for(int i = 0; i < partitionThreadCount; ++i)
        {
            PartitionThread partitionThread = new PartitionThread(mConfig, inputBamFiles, taskQueue, fileWriterCache);
            partitionThreads.add(partitionThread);
        }

        return partitionThreads;
    }

    private void setReadLength()
    {
        if(mConfig.readLength() > 0) // skip if set in config
            return;

        // sample the BAM to determine read length
        BamSampler bamSampler = new BamSampler(mConfig.RefGenomeFile);

        ChrBaseRegion sampleRegion = !mConfig.SpecificChrRegions.Regions.isEmpty() ?
                mConfig.SpecificChrRegions.Regions.get(0) : bamSampler.defaultRegion();

        int readLength = DEFAULT_READ_LENGTH;

        if(bamSampler.calcBamCharacteristics(mConfig.BamFiles.get(0), sampleRegion) && bamSampler.maxReadLength() > 0)
        {
            readLength = bamSampler.maxReadLength();
            RD_LOGGER.debug("BAM sampled max read-length({})", readLength);
        }
        else
        {
            RD_LOGGER.debug("BAM read-length sampling failed, using default read length({})", DEFAULT_READ_LENGTH);
        }

        mConfig.setReadLength(readLength);
    }

    private List<PerformanceCounter> mergePerfCounters(final List<PartitionReader> partitionReaders)
    {
        List<PerformanceCounter> combinedPerfCounters = partitionReaders.get(0).perfCounters();

        for(int i = 1; i < partitionReaders.size(); ++i)
        {
            List<PerformanceCounter> chrPerfCounters = partitionReaders.get(i).perfCounters();

            for(int j = 0; j < chrPerfCounters.size(); ++j)
            {
                combinedPerfCounters.get(j).merge(chrPerfCounters.get(j));
            }
        }

        return combinedPerfCounters;
    }

    private void logPerformanceStats(final List<PerformanceCounter> combinedPerfCounters)
    {
        if(mConfig.perfDebug())
        {
            for(int j = 0; j < combinedPerfCounters.size(); ++j)
            {
                PerformanceCounter perfCounter = combinedPerfCounters.get(j);

                if(j == 0)
                    perfCounter.logIntervalStats(10);
                else
                    perfCounter.logStats();
            }
        }
        else
        {
            combinedPerfCounters.forEach(x -> x.logStats());
        }
    }

    public static void main(final String[] args)
    {
        ConfigBuilder configBuilder = new ConfigBuilder(APP_NAME);
        registerConfig(configBuilder);

        configBuilder.checkAndParseCommandLine(args);

        ReduxApplication reduxApplication = new ReduxApplication(configBuilder);
        reduxApplication.run();
    }
}
