package com.hartwig.hmftools.redux;

import static java.lang.Math.ceil;
import static java.lang.Math.round;

import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeCoordinates.refGenomeCoordinates;
import static com.hartwig.hmftools.common.utils.PerformanceCounter.runTimeMinsStr;
import static com.hartwig.hmftools.redux.ReduxConfig.APP_NAME;
import static com.hartwig.hmftools.redux.ReduxConfig.RD_LOGGER;
import static com.hartwig.hmftools.redux.ReduxConfig.registerConfig;
import static com.hartwig.hmftools.redux.common.Constants.DEFAULT_READ_LENGTH;
import static com.hartwig.hmftools.redux.common.ReadUnmapper.unmapMateAlignment;
import static com.hartwig.hmftools.redux.common.ReadUnmapper.unmapReadAlignment;

import java.util.Collections;
import java.util.List;
import java.util.concurrent.Callable;
import java.util.concurrent.atomic.AtomicLong;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.bamops.BamSampler;
import com.hartwig.hmftools.common.basequal.jitter.ConsensusMarker;
import com.hartwig.hmftools.common.basequal.jitter.JitterAnalyser;
import com.hartwig.hmftools.common.basequal.jitter.JitterAnalyserConfig;
import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeCoordinates;
import com.hartwig.hmftools.common.region.ChrBaseRegion;
import com.hartwig.hmftools.common.utils.PerformanceCounter;
import com.hartwig.hmftools.common.utils.TaskExecutor;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;
import com.hartwig.hmftools.redux.common.FragmentStatus;
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
                    mConfig.JitterMsiMaxSitePercContribution);

            ConsensusMarker consensusMarker = ConsensusMarker.fromSequencingType(mConfig.Sequencing);
            jitterAnalyser = new JitterAnalyser(jitterConfig, RD_LOGGER, consensusMarker);
        }

        FileWriterCache fileWriterCache = new FileWriterCache(mConfig, jitterAnalyser);

        // partition the genome into sequential regions to be processed by each thread
        List<PartitionReader> partitionReaders = allocateGenomeRegions(fileWriterCache);

        List<Callable> callableTasks = partitionReaders.stream().collect(Collectors.toList());

        RD_LOGGER.debug("splitting regions across {} threads", partitionReaders.size());
        if(!TaskExecutor.executeTasks(callableTasks, mConfig.Threads))
            System.exit(1);

        RD_LOGGER.info("all partition tasks complete");

        /*
        // write any orphaned or remaining fragments (can be supplementaries)
        BamWriter recordWriter = fileWriterCache.getUnsortedBamWriter();
        */

        ConsensusReads consensusReads = new ConsensusReads(mConfig.RefGenome);
        consensusReads.setDebugOptions(mConfig.RunChecks);

        long unmappedReads = writeUnmappedReads(fileWriterCache);

        int totalUnwrittenFragments = 0;

        if(totalUnwrittenFragments > 0)
        {
            RD_LOGGER.info("wrote {} remaining cached fragments", totalUnwrittenFragments);
        }

        List<PerformanceCounter> combinedPerfCounters = mergePerfCounters(partitionReaders);

        Statistics combinedStats = new Statistics();
        partitionReaders.forEach(x -> combinedStats.merge(x.statistics()));
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

        if(!mConfig.JitterMsiOnly)
        {
            combinedStats.logStats();

            long totalWrittenReads = fileWriterCache.totalWrittenReads();
            long unmappedDroppedReads = mConfig.UnmapRegions.stats().SupplementaryCount.get() + mConfig.UnmapRegions.stats().SecondaryCount.get();

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
        }

        logPerformanceStats(combinedPerfCounters);

        RD_LOGGER.info("Redux complete, mins({})", runTimeMinsStr(startTimeMs));
    }

    private List<PartitionReader> allocateGenomeRegions(final FileWriterCache fileWriterCache)
    {
        List<List<ChrBaseRegion>> partitionRegions = Lists.newArrayList();

        if(!mConfig.SpecificChrRegions.Regions.isEmpty())
        {
            // could split long regions but for now just split by raw count across threads
            int regionCount = mConfig.SpecificChrRegions.Regions.size();
            int regionsPerThread = (int)round(regionCount / mConfig.Threads);

            List<ChrBaseRegion> currentRegions = null;
            for(ChrBaseRegion region : mConfig.SpecificChrRegions.Regions)
            {
                if(currentRegions == null || currentRegions.size() >= regionsPerThread)
                {
                    currentRegions = Lists.newArrayList();
                    partitionRegions.add(currentRegions);
                }

                currentRegions.add(region);
            }
        }
        else
        {
            // divide the genome into equal lengths
            RefGenomeCoordinates refGenomeCoordinates = refGenomeCoordinates(mConfig.RefGenVersion);
            long totalLength = 0;

            for(HumanChromosome chromosome : HumanChromosome.values())
            {
                String chromosomeStr = mConfig.RefGenVersion.versionedChromosome(chromosome.toString());

                if(mConfig.SpecificChrRegions.excludeChromosome(chromosomeStr))
                    continue;

                totalLength += refGenomeCoordinates.Lengths.get(chromosome);
            }

            long intervalLength = (int)ceil(totalLength / (double)mConfig.Threads);
            int chrEndBuffer = (int)round(intervalLength * 0.05);
            long currentLength = 0;
            int nextRegionStart = 1;
            List<ChrBaseRegion> currentRegions = Lists.newArrayList();
            partitionRegions.add(currentRegions);

            for(HumanChromosome chromosome : HumanChromosome.values())
            {
                String chromosomeStr = mConfig.RefGenVersion.versionedChromosome(chromosome.toString());

                if(mConfig.SpecificChrRegions.excludeChromosome(chromosomeStr))
                    continue;

                nextRegionStart = 1;
                int chromosomeLength = refGenomeCoordinates.Lengths.get(chromosome);
                int remainingChromosomeLength = chromosomeLength;

                while(currentLength + remainingChromosomeLength >= intervalLength)
                {
                    int remainingIntervalLength = (int)(intervalLength - currentLength);
                    int regionEnd = nextRegionStart + remainingIntervalLength - 1;

                    if(chromosomeLength - regionEnd < chrEndBuffer)
                    {
                        regionEnd = chromosomeLength;
                        remainingChromosomeLength = 0;
                    }

                    currentRegions.add(new ChrBaseRegion(chromosomeStr, nextRegionStart, regionEnd));

                    currentRegions = Lists.newArrayList();
                    partitionRegions.add(currentRegions);
                    currentLength = 0;

                    if(remainingChromosomeLength == 0)
                        break;

                    nextRegionStart = regionEnd + 1;
                    remainingChromosomeLength = chromosomeLength - nextRegionStart + 1;
                }

                if(remainingChromosomeLength <= 0)
                    continue;

                currentLength += remainingChromosomeLength;

                currentRegions.add(new ChrBaseRegion(chromosomeStr, nextRegionStart, chromosomeLength));
            }
        }

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
            PartitionReader partitionReader = new PartitionReader(mConfig, regions, bamWriter);
            partitionReaders.add(partitionReader);
        }

        return partitionReaders;
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

    private void logPerformanceStats(final List<PerformanceCounter> combinedPerfCounters)
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
