package com.hartwig.hmftools.redux;

import static java.lang.Math.max;
import static java.lang.Math.min;
import static java.lang.String.format;

import static com.hartwig.hmftools.common.region.PartitionUtils.partitionChromosome;
import static com.hartwig.hmftools.common.utils.PerformanceCounter.runTimeMinsStr;
import static com.hartwig.hmftools.common.utils.TaskExecutor.runThreadTasks;
import static com.hartwig.hmftools.redux.ReduxConfig.APP_NAME;
import static com.hartwig.hmftools.redux.ReduxConfig.RD_LOGGER;
import static com.hartwig.hmftools.redux.ReduxConfig.registerConfig;
import static com.hartwig.hmftools.redux.common.Constants.DEFAULT_READ_LENGTH;
import static com.hartwig.hmftools.redux.common.Constants.LOCK_ACQUIRE_LONG_TIME_MS;
import static com.hartwig.hmftools.redux.common.ReadUnmapper.unmapMateAlignment;
import static com.hartwig.hmftools.redux.common.ReadUnmapper.unmapReadAlignment;

import java.util.ArrayList;
import java.util.List;
import java.util.Queue;
import java.util.concurrent.ConcurrentLinkedQueue;
import java.util.concurrent.atomic.AtomicLong;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.bam.BamSampler;
import com.hartwig.hmftools.common.basequal.jitter.JitterAnalyser;
import com.hartwig.hmftools.common.basequal.jitter.JitterAnalyserConfig;
import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.region.ChrBaseRegion;
import com.hartwig.hmftools.common.utils.PerformanceCounter;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;
import com.hartwig.hmftools.redux.common.FragmentStatus;
import com.hartwig.hmftools.redux.common.PartitionData;
import com.hartwig.hmftools.redux.common.Statistics;
import com.hartwig.hmftools.redux.consensus.ConsensusReads;
import com.hartwig.hmftools.redux.write.BamWriter;
import com.hartwig.hmftools.redux.write.FileWriterCache;

import htsjdk.samtools.SAMRecord;

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
                    JitterAnalyserConfig.DEFAULT_MIN_MAPPING_QUALITY, mConfig.JitterMaxSitesPerType, mConfig.Threads, false);

            jitterAnalyser = new JitterAnalyser(jitterConfig, RD_LOGGER);
        }

        FileWriterCache fileWriterCache = new FileWriterCache(mConfig, jitterAnalyser);

        PartitionDataStore partitionDataStore = new PartitionDataStore(mConfig);

        // partition all chromosomes
        Queue<ChrBaseRegion> partitions = new ConcurrentLinkedQueue<>();

        for(HumanChromosome chromosome : HumanChromosome.values())
        {
            String chromosomeStr = mConfig.RefGenVersion.versionedChromosome(chromosome.toString());

            if(mConfig.SpecificChrRegions.excludeChromosome(chromosomeStr))
                continue;

            List<ChrBaseRegion> chrPartitions = partitionChromosome(
                    chromosomeStr, mConfig.RefGenVersion, mConfig.SpecificChrRegions.Regions, mConfig.PartitionSize);

            partitions.addAll(chrPartitions);
        }

        List<PartitionThread> partitionTasks = Lists.newArrayList();
        List<Thread> workers = new ArrayList<>();

        int partitionCount = partitions.size();

        for(int i = 0; i < min(partitionCount, mConfig.Threads); ++i)
        {
            PartitionThread partitionThread = new PartitionThread(i, mConfig, partitions, fileWriterCache, partitionDataStore);
            partitionTasks.add(partitionThread);
            workers.add(partitionThread);
        }

        RD_LOGGER.debug("splitting {} partitions across {} threads", partitionCount, partitionTasks.size());

        if(!runThreadTasks(workers))
            System.exit(1);

        RD_LOGGER.info("all partition tasks complete");

        long unmappedReads = writeUnmappedReads(fileWriterCache);

        List<PartitionReader> partitionReaders = partitionTasks.stream().map(x -> x.partitionReader()).collect(Collectors.toList());

        int maxLogFragments = (mConfig.RunChecks || mConfig.LogFinalCache) ? 100 : 0;
        int totalUnwrittenFragments = 0;
        ConsensusReads consensusReads = new ConsensusReads(mConfig.RefGenome);
        consensusReads.setDebugOptions(mConfig.RunChecks);

        // write any orphaned or remaining fragments (can be supplementaries)
        BamWriter recordWriter = fileWriterCache.getUnsortedBamWriter();

        for(PartitionData partitionData : partitionDataStore.partitions())
        {
            int cachedReadCount = partitionData.writeRemainingReads(recordWriter, consensusReads, maxLogFragments > 0);
            totalUnwrittenFragments += cachedReadCount;
            maxLogFragments = max(0, maxLogFragments - cachedReadCount);
        }

        if(totalUnwrittenFragments > 0)
        {
            RD_LOGGER.info("wrote {} remaining cached fragments", totalUnwrittenFragments);
        }

        List<PerformanceCounter> combinedPerfCounters = mergePerfCounters(partitionReaders);

        Statistics combinedStats = new Statistics();
        partitionReaders.forEach(x -> combinedStats.merge(x.statistics()));
        partitionDataStore.partitions().forEach(x -> combinedStats.merge(x.statistics()));
        combinedStats.ConsensusStats.merge(consensusReads.consensusStats());

        // free up any processing state
        partitionReaders.clear();

        fileWriterCache.close();

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

        if(fileWriterCache.runSortMergeIndex())
        {
            // usually avoid manual calls to this but since the external BAM tools make independent calls to access memory and
            // the core routines are complete, it is helpful to do so now
            System.gc();

            // log interim time
            RD_LOGGER.info("BAM duplicate processing complete, mins({})", runTimeMinsStr(startTimeMs));

            if(!fileWriterCache.sortAndIndexBams())
            {
                RD_LOGGER.error("sort-merge-index failed");
                System.exit(1);
            }
        }

        combinedStats.logStats();

        long totalWrittenReads = fileWriterCache.totalWrittenReads();
        long unmappedDroppedReads = mConfig.UnmapRegions.stats().SupplementaryCount.get();

        if(mConfig.UnmapRegions.enabled())
        {
            RD_LOGGER.info("unmapped stats: {}", mConfig.UnmapRegions.stats().toString());
        }

        if(combinedStats.TotalReads + unmappedReads != totalWrittenReads + unmappedDroppedReads)
        {
            long difference = combinedStats.TotalReads + unmappedReads - totalWrittenReads - unmappedDroppedReads;

            RD_LOGGER.warn("reads processed({}) vs written({}) mismatch diffLessDropped({})",
                    combinedStats.TotalReads + unmappedReads, totalWrittenReads, difference);
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

        logPerformanceStats(combinedPerfCounters, partitionDataStore);

        RD_LOGGER.info("Mark duplicates complete, mins({})", runTimeMinsStr(startTimeMs));
    }

    private long writeUnmappedReads(final FileWriterCache fileWriterCache)
    {
        if(mConfig.SpecificChrRegions.hasFilters() || !mConfig.WriteBam)
            return 0;

        BamWriter bamWriter = fileWriterCache.getUnsortedBamWriter();

        BamReader bamReader = new BamReader(mConfig);

        AtomicLong unmappedCount = new AtomicLong();
        AtomicLong nonHumanContigCount = new AtomicLong();

        bamReader.queryUnmappedReads((final SAMRecord record) ->
        {
            bamWriter.writeRead(record, FragmentStatus.UNSET);
            unmappedCount.incrementAndGet();
        });

        // do the same for non-human contigs
        bamReader.queryNonHumanContigs((final SAMRecord record) ->
        {
            processNonHumanContigReads(record, bamWriter);
            nonHumanContigCount.incrementAndGet();
        });

        if(unmappedCount.get() > 0 || nonHumanContigCount.get() > 0)
        {
            RD_LOGGER.debug("wrote unmapped({}) otherContig({}) reads", unmappedCount, nonHumanContigCount);
        }

        return unmappedCount.get() + nonHumanContigCount.get();
    }

    private void processNonHumanContigReads(final SAMRecord record, final BamWriter bamWriter)
    {
        // if these have a mate in a human chromosome, then they have been unmapped in that read, so do so here as well
        if(record.getReadPairedFlag() && !record.getMateUnmappedFlag() && HumanChromosome.contains(record.getMateReferenceName()))
        {
            if(record.getSupplementaryAlignmentFlag())
                return; // drop as per standard logic

            boolean mateUnmapped = mConfig.UnmapRegions.mateInUnmapRegion(record);

            // if the human-chromosome mate was unmapped (ie in an unmap region), then this read should also now be unmapped
            // otherwise it should be unmapped but leave its mate attributes as-is
            unmapReadAlignment(record, mateUnmapped, mateUnmapped);

            if(mateUnmapped)
            {
                unmapMateAlignment(record, false, true);
            }
        }

        bamWriter.writeRead(record, FragmentStatus.UNSET);
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

    private void logPerformanceStats(final List<PerformanceCounter> combinedPerfCounters, final PartitionDataStore partitionDataStore)
    {
        if(mConfig.PerfDebug)
        {
            for(int j = 0; j < combinedPerfCounters.size(); ++j)
            {
                PerformanceCounter perfCounter = combinedPerfCounters.get(j);

                if(j == 0)
                    perfCounter.logIntervalStats(10);
                else
                    perfCounter.logStats();
            }

            // check partition store locking times
            double totalLockTimeMs = 0;

            for(PartitionData partitionData : partitionDataStore.partitions())
            {
                double lockTime = partitionData.totalLockTimeMs();

                totalLockTimeMs += lockTime;

                if(lockTime > LOCK_ACQUIRE_LONG_TIME_MS)
                {
                    RD_LOGGER.debug("partition({}) lock-acquisition time({}ms)",
                            partitionData.partitionStr(), format("%.1f", lockTime));
                }
            }

            if(totalLockTimeMs > LOCK_ACQUIRE_LONG_TIME_MS)
            {
                RD_LOGGER.debug("partition cache total lock-acquisition time({}s)",
                        format("%.3f", totalLockTimeMs / 1000));
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
