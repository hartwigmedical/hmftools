package com.hartwig.hmftools.markdups;

import static java.lang.Math.max;
import static java.lang.String.format;

import static com.hartwig.hmftools.common.utils.ConfigUtils.setLogLevel;
import static com.hartwig.hmftools.markdups.MarkDupsConfig.MD_LOGGER;

import java.util.Collections;
import java.util.List;
import java.util.concurrent.Callable;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeCoordinates;
import com.hartwig.hmftools.common.utils.PerformanceCounter;
import com.hartwig.hmftools.common.utils.TaskExecutor;
import com.hartwig.hmftools.common.utils.sv.ChrBaseRegion;
import com.hartwig.hmftools.common.utils.version.VersionInfo;
import com.hartwig.hmftools.markdups.common.Fragment;
import com.hartwig.hmftools.markdups.common.PartitionData;
import com.hartwig.hmftools.markdups.common.Statistics;
import com.hartwig.hmftools.markdups.umi.ConsensusReads;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.jetbrains.annotations.NotNull;

public class MarkDuplicates
{
    private final MarkDupsConfig mConfig;

    public MarkDuplicates(final CommandLine cmd)
    {
        mConfig = new MarkDupsConfig(cmd);
    }

    public void run()
    {
        if(!mConfig.isValid())
            System.exit(1);

        MD_LOGGER.info("sample({}) starting mark duplicates", mConfig.SampleId);

        long startTimeMs = System.currentTimeMillis();

        List<ChromosomeReader> chromosomeReaders = Lists.newArrayList();

        RefGenomeCoordinates refGenomeCoordinates = mConfig.RefGenVersion.is37() ? RefGenomeCoordinates.COORDS_37 : RefGenomeCoordinates.COORDS_38;

        RecordWriter recordWriter = new RecordWriter(mConfig);
        PartitionDataStore partitionDataStore = new PartitionDataStore(mConfig);
        final List<Callable> callableList = Lists.newArrayList();

        for(HumanChromosome chromosome : HumanChromosome.values())
        {
            String chromosomeStr = mConfig.RefGenVersion.versionedChromosome(chromosome.toString());

            if(!mConfig.SpecificChromosomes.isEmpty() && !mConfig.SpecificChromosomes.contains(chromosomeStr))
                continue;

            ChrBaseRegion chrBaseRegion = new ChrBaseRegion(chromosomeStr, 1, refGenomeCoordinates.Lengths.get(chromosome));

            ChromosomeReader chromosomeReader = new ChromosomeReader(chrBaseRegion, mConfig, recordWriter, partitionDataStore);
            chromosomeReaders.add(chromosomeReader);
            callableList.add(chromosomeReader);
        }

        if(!TaskExecutor.executeTasks(callableList, mConfig.Threads))
            System.exit(1);

        int maxLogFragments = (mConfig.RunChecks || mConfig.PerfDebug) ? 100 : 0;
        int totalUnwrittenFragments = 0;
        ConsensusReads consensusReads = new ConsensusReads(mConfig.UMIs, mConfig.RefGenome);
        for(PartitionData partitionData : partitionDataStore.partitions())
        {
            int cachedReadCount = partitionData.writeRemainingReads(recordWriter, consensusReads, maxLogFragments > 0);
            totalUnwrittenFragments += cachedReadCount;
            maxLogFragments = max(0, maxLogFragments - cachedReadCount);
        }

        if(totalUnwrittenFragments > 0)
        {
            MD_LOGGER.info("wrote {} cached fragments", totalUnwrittenFragments);
        }

        recordWriter.close();

        MD_LOGGER.debug("all chromosome tasks complete");

        Statistics combinedStats = new Statistics();
        chromosomeReaders.forEach(x -> combinedStats.merge(x.statistics()));
        partitionDataStore.partitions().forEach(x -> combinedStats.merge(x.statistics()));

        combinedStats.logStats();

        if(combinedStats.TotalReads != recordWriter.recordWriteCount())
        {
            MD_LOGGER.warn("reads processed({}) vs written({}) mismatch", combinedStats.TotalReads, recordWriter.recordWriteCount());
            recordWriter.logUnwrittenReads();
        }

        if(mConfig.UMIs.Debug)
        {
            combinedStats.writeUmiStats(mConfig);
        }

        PerformanceCounter combinedPerfCounter = chromosomeReaders.get(0).perfCounter();

        for(int i = 1; i < chromosomeReaders.size(); ++i)
        {
            combinedPerfCounter.merge(chromosomeReaders.get(i).perfCounter());
        }

        if(mConfig.PerfDebug)
        {
            combinedPerfCounter.logIntervalStats(10);

            List<Double> partitionLockTimes = Lists.newArrayList();
            partitionDataStore.partitions().forEach(x -> partitionLockTimes.add(x.totalLockTime()));
            Collections.sort(partitionLockTimes, Collections.reverseOrder());
            double nthTime = partitionLockTimes.size() >= 5 ? partitionLockTimes.get(4) : partitionLockTimes.get(partitionLockTimes.size() - 1);

            for(PartitionData partitionData : partitionDataStore.partitions())
            {
                double lockTime = partitionData.totalLockTime();

                if(lockTime > 0 && lockTime >= nthTime)
                {
                    MD_LOGGER.debug("partition({}) total lock-acquisition time({})",
                            partitionData.partitionStr(), format("%.3f", lockTime));
                }
            }
        }
        else
        {
            combinedPerfCounter.logStats();
        }

        long timeTakenMs = System.currentTimeMillis() - startTimeMs;
        double timeTakeMins = timeTakenMs / 60000.0;

        MD_LOGGER.info("Mark duplicates complete, mins({})", format("%.3f", timeTakeMins));
    }

    public static void main(@NotNull final String[] args)
    {
        final VersionInfo version = new VersionInfo("mark-dups.version");
        MD_LOGGER.info("MarkDups version: {}", version.version());

        final Options options = MarkDupsConfig.createCmdLineOptions();

        try
        {
            final CommandLine cmd = createCommandLine(args, options);

            setLogLevel(cmd);

            MarkDuplicates markDuplicates = new MarkDuplicates(cmd);
            markDuplicates.run();
        }
        catch(ParseException e)
        {
            MD_LOGGER.warn(e);
            final HelpFormatter formatter = new HelpFormatter();
            formatter.printHelp("MarkDuplicates", options);
            System.exit(1);
        }
    }

    @NotNull
    private static CommandLine createCommandLine(@NotNull final String[] args, @NotNull final Options options) throws ParseException
    {
        final CommandLineParser parser = new DefaultParser();
        return parser.parse(options, args);
    }
}
