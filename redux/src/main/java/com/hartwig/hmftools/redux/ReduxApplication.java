package com.hartwig.hmftools.redux;

import static com.hartwig.hmftools.common.utils.PerformanceCounter.runTimeMinsStr;
import static com.hartwig.hmftools.redux.ReduxConfig.APP_NAME;
import static com.hartwig.hmftools.redux.ReduxConfig.RD_LOGGER;
import static com.hartwig.hmftools.redux.ReduxConfig.registerConfig;
import static com.hartwig.hmftools.redux.ReduxConfig.splitRegionsByThreads;
import static com.hartwig.hmftools.redux.common.Constants.DEFAULT_READ_LENGTH;

import java.util.Collections;
import java.util.List;
import java.util.concurrent.Callable;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.bamops.BamSampler;
import com.hartwig.hmftools.common.basequal.jitter.JitterAnalyser;
import com.hartwig.hmftools.common.basequal.jitter.JitterAnalyserConfig;
import com.hartwig.hmftools.common.region.ChrBaseRegion;
import com.hartwig.hmftools.common.utils.PerformanceCounter;
import com.hartwig.hmftools.common.utils.TaskExecutor;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;
import com.hartwig.hmftools.redux.common.Statistics;
import com.hartwig.hmftools.redux.write.BamWriter;
import com.hartwig.hmftools.redux.write.FileWriterCache;
import com.hartwig.hmftools.redux.write.FinalBamProcessor;

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
                    mConfig.JitterMsiMaxSitePercContribution);

            jitterAnalyser = new JitterAnalyser(jitterConfig, RD_LOGGER);
        }

        FileWriterCache fileWriterCache = new FileWriterCache(mConfig, jitterAnalyser);

        // partition the genome into sequential regions to be processed by each thread
        List<PartitionReader> partitionReaders = allocateGenomeRegions(fileWriterCache);

        List<Callable> callableTasks = partitionReaders.stream().collect(Collectors.toList());

        RD_LOGGER.debug("splitting regions across {} threads", partitionReaders.size());
        if(!TaskExecutor.executeTasks(callableTasks, mConfig.Threads))
            System.exit(1);

        RD_LOGGER.info("all partition tasks complete");

        // TODO: consider leaving the read TSV writer open so it can write unmapped read info
        fileWriterCache.close();

        Statistics combinedStats = new Statistics();
        partitionReaders.forEach(x -> combinedStats.merge(x.statistics()));

        if(mConfig.BamToolPath != null)
        {
            // usually avoid manual calls to this but since the external BAM tools make independent calls to access memory and
            // the core routines are complete, it is helpful to do so now
            System.gc();

            // log interim time
            RD_LOGGER.info("BAM duplicate processing complete, mins({})", runTimeMinsStr(startTimeMs));

            FinalBamProcessor finalBamProcessor = new FinalBamProcessor(mConfig, fileWriterCache);
            if(!finalBamProcessor.run())
                System.exit(1);

            combinedStats.ConsensusStats.merge(finalBamProcessor.statistics().ConsensusStats);
        }

        long unmappedReads = 0; // writeUnmappedReads(fileWriterCache);

        List<PerformanceCounter> combinedPerfCounters = mergePerfCounters(partitionReaders);


        // free up any processing state
        partitionReaders.clear();


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

        if(!mConfig.JitterMsiOnly)
        {
            combinedStats.logStats();

            long totalWrittenReads = fileWriterCache.totalWrittenReads();
            long unmappedDroppedReads = mConfig.UnmapRegions.stats().SupplementaryCount.get() + mConfig.UnmapRegions.stats().SecondaryCount.get();

            if(mConfig.UnmapRegions.enabled())
            {
                RD_LOGGER.info("unmapped stats: {}", mConfig.UnmapRegions.stats().toString());
            }

            /*
            if(combinedStats.TotalReads + unmappedReads != totalWrittenReads + unmappedDroppedReads)
            {
                long difference = combinedStats.TotalReads + unmappedReads - totalWrittenReads - unmappedDroppedReads;

                RD_LOGGER.warn("reads processed({}) vs written({}) mismatch diffLessDropped({})",
                        combinedStats.TotalReads + unmappedReads, totalWrittenReads, difference);
            }
            */

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
        }

        logPerformanceStats(combinedPerfCounters);

        RD_LOGGER.info("Redux complete, mins({})", runTimeMinsStr(startTimeMs));
    }

    private List<PartitionReader> allocateGenomeRegions(final FileWriterCache fileWriterCache)
    {
        List<List<ChrBaseRegion>> partitionRegions = splitRegionsByThreads(mConfig.SpecificChrRegions, mConfig.Threads, mConfig.RefGenVersion);

        if(partitionRegions.isEmpty())
            return Collections.emptyList();

        List<PartitionReader> partitionReaders = Lists.newArrayListWithCapacity(mConfig.Threads);
        int threadIndex = 0;

        for(List<ChrBaseRegion> regions : partitionRegions)
        {
            if(regions.isEmpty())
                break;

            long regionsLength = regions.stream().mapToLong(x -> x.baseLength()).sum();
            String regionsStr = regions.stream().map(x -> x.toString()).collect(Collectors.joining(";"));
            RD_LOGGER.debug("adding partition regions({}) totalLength({}): {}}", regions.size(), regionsLength, regionsStr);

            BamWriter bamWriter = fileWriterCache.getPartitionBamWriter(String.valueOf(threadIndex++));

            PartitionReader partitionReader = new PartitionReader(
                    mConfig, regions, mConfig.BamFiles, bamWriter, fileWriterCache.getUnsortedBamWriter());

            partitionReaders.add(partitionReader);
        }

        return partitionReaders;
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
