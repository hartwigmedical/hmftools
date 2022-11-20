package com.hartwig.hmftools.bamtools;

import static java.lang.Math.min;
import static java.lang.String.format;

import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeFunctions.stripChrPrefix;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion.V37;
import static com.hartwig.hmftools.common.utils.ConfigUtils.setLogLevel;

import java.util.ArrayList;
import java.util.List;
import java.util.Queue;
import java.util.concurrent.ConcurrentLinkedQueue;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.bamtools.metrics.CombinedStats;
import com.hartwig.hmftools.bamtools.metrics.MetricsWriter;
import com.hartwig.hmftools.bamtools.slice.SliceWriter;
import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeCoordinates;
import com.hartwig.hmftools.common.utils.sv.ChrBaseRegion;
import com.hartwig.hmftools.common.utils.version.VersionInfo;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.jetbrains.annotations.NotNull;

public class BamToolsApplication
{
    private final BmConfig mConfig;

    public BamToolsApplication(final CommandLine cmd)
    {
        mConfig = new BmConfig(cmd);
    }

    public void run()
    {
        if(!mConfig.isValid())
            System.exit(1);

        BmConfig.BM_LOGGER.info("sample({}) starting bam metrics", mConfig.SampleId);

        long startTimeMs = System.currentTimeMillis();

        List<ChrBaseRegion> allRegions = Lists.newArrayList();

        for(HumanChromosome chromosome : HumanChromosome.values())
        {
            String chromosomeStr = mConfig.RefGenVersion.versionedChromosome(chromosome.toString());

            if(!mConfig.SpecificChromosomes.isEmpty() && !mConfig.SpecificChromosomes.contains(chromosomeStr))
                continue;

            allRegions.addAll(partitionChromosome(chromosomeStr));
        }

        BmConfig.BM_LOGGER.info("splitting {} regions across {} threads", allRegions.size(), mConfig.Threads);

        Queue<PartitionTask> partitions = new ConcurrentLinkedQueue<>();

        int taskId = 0;
        for(int i = 0; i < allRegions.size(); ++i)
        {
            partitions.add(new PartitionTask(allRegions.get(i), taskId++));
        }

        CombinedStats combinedStats = new CombinedStats(mConfig.MaxCoverage);

        SliceWriter sliceWriter = new SliceWriter(mConfig);

        List<Thread> workers = new ArrayList<>();

        for(int i = 0; i < min(allRegions.size(), mConfig.Threads); ++i)
        {
            workers.add(new PartitionThread(mConfig, partitions, combinedStats, sliceWriter));
        }

        for(Thread worker : workers)
        {
            try
            {
                worker.join();
            }
            catch(InterruptedException e)
            {
                BmConfig.BM_LOGGER.error("task execution error: {}", e.toString());
                e.printStackTrace();
            }
        }

        combinedStats.metrics().finalise(mConfig.ExcludeZeroCoverage);
        sliceWriter.close();

        BmConfig.BM_LOGGER.info("all regions complete, totalReads({}) stats: {}", combinedStats.totalReads(), combinedStats.metrics());

        if(mConfig.PerfDebug)
            combinedStats.perfCounter().logIntervalStats(10);
        else
            combinedStats.perfCounter().logStats();

        MetricsWriter.writeResults(combinedStats.metrics(), mConfig);

        long timeTakenMs = System.currentTimeMillis() - startTimeMs;
        double timeTakeMins = timeTakenMs / 60000.0;

        BmConfig.BM_LOGGER.info("BamMetrics complete, mins({})", format("%.3f", timeTakeMins));
    }

    private List<ChrBaseRegion> partitionChromosome(final String chromosome)
    {
        if(!mConfig.SpecificRegions.isEmpty())
        {
            List<ChrBaseRegion> partitions = Lists.newArrayList();

            for(ChrBaseRegion region : mConfig.SpecificRegions)
            {
                if(region.Chromosome.equals(chromosome))
                {
                    partitions.addAll(buildPartitions(chromosome, region.start() ,region.end()));
                }
            }

            return partitions;
        }

        RefGenomeCoordinates refGenomeCoords = mConfig.RefGenVersion == V37 ? RefGenomeCoordinates.COORDS_37 : RefGenomeCoordinates.COORDS_38;
        int chromosomeLength = refGenomeCoords.length(stripChrPrefix(chromosome));
        return buildPartitions(chromosome, 1, chromosomeLength);
    }

    private List<ChrBaseRegion> buildPartitions(final String chromosome, int minPosition, int maxPosition)
    {
        final List<ChrBaseRegion> partitions = Lists.newArrayList();

        for(int i = 0; ; i++)
        {
            int start = minPosition + i * mConfig.PartitionSize;
            int end = min(start + mConfig.PartitionSize - 1, maxPosition);
            partitions.add(new ChrBaseRegion(chromosome, start, end));

            if(end >= maxPosition)
                break;
        }

        return partitions;
    }


    public static void main(@NotNull final String[] args) throws Exception
    {
        final VersionInfo version = new VersionInfo("bam-tools.version");
        BmConfig.BM_LOGGER.info("BamMetrics version: {}", version.version());

        final Options options = BmConfig.createCmdLineOptions();

        try
        {
            final CommandLine cmd = createCommandLine(args, options);

            setLogLevel(cmd);

            BamToolsApplication bamMetrics = new BamToolsApplication(cmd);
            bamMetrics.run();
        }
        catch(ParseException e)
        {
            BmConfig.BM_LOGGER.warn(e);
            final HelpFormatter formatter = new HelpFormatter();
            formatter.printHelp("BamMetrics", options);
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
