package com.hartwig.hmftools.bammetrics;

import static java.lang.String.format;

import static com.hartwig.hmftools.bammetrics.BmConfig.BM_LOGGER;
import static com.hartwig.hmftools.bammetrics.BmConfig.createCmdLineOptions;
import static com.hartwig.hmftools.common.utils.ConfigUtils.setLogLevel;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.createBufferedWriter;

import java.io.BufferedWriter;
import java.io.IOException;
import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.utils.version.VersionInfo;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.jetbrains.annotations.NotNull;

public class BamMetricsApplication
{
    private final BmConfig mConfig;

    public BamMetricsApplication(final CommandLine cmd)
    {
        mConfig = new BmConfig(cmd);
    }

    public void run()
    {
        if(!mConfig.isValid())
            System.exit(1);

        BM_LOGGER.info("sample({}) starting bam metrics", mConfig.SampleId);

        long startTimeMs = System.currentTimeMillis();

        Metrics combinedMetrics = new Metrics(mConfig.MaxCoverage);

        for(HumanChromosome chromosome : HumanChromosome.values())
        {
            String chromosomeStr = mConfig.RefGenVersion.versionedChromosome(chromosome.toString());

            if(!mConfig.SpecificChromosomes.isEmpty() && !mConfig.SpecificChromosomes.contains(chromosomeStr))
                continue;

            BM_LOGGER.info("processing chromosome({})", chromosomeStr);

            ChromosomeTask chromosomeTask = new ChromosomeTask(chromosomeStr, mConfig);
            chromosomeTask.process();
            combinedMetrics.merge(chromosomeTask.combinedMetrics());

            if(mConfig.PerfDebug)
            {
                // combinedStats.PerfCounters.get(i).merge(chromosomeTask.combinedMetrics().PerfCounters.get((i)));
            }

            System.gc();
        }

        combinedMetrics.finalise(mConfig.ExcludeZeroCoverage);
        writeResults(combinedMetrics);

        long timeTakenMs = System.currentTimeMillis() - startTimeMs;
        double timeTakeMins = timeTakenMs / 60000.0;

        if(mConfig.PerfDebug)
        {
            // combinedStats.PerfCounters.forEach(x -> x.logStats());
        }

        BM_LOGGER.info("BamMetrics complete, mins({})", format("%.3f", timeTakeMins));
    }

    private void writeResults(final Metrics metrics)
    {
        try
        {
            // write summary metrics
            String filename = mConfig.formFilename("metrics");
            BufferedWriter writer = createBufferedWriter(filename, false);

            writer.write("GenomeTerritory,CoverageTotal,CoverageMean,CoverageMedian,CoverageStdDev");

            for(FilterType type : FilterType.values())
            {
                if(type == FilterType.UNFILTERED)
                    continue;

                writer.write(format(",Pct%s", type.description()));
            }

            List<Integer> coverageLevels = Lists.newArrayList(1, 5, 10, 15, 20, 25, 30, 40, 50, 60, 70, 80, 90, 100);

            for(Integer coverage : coverageLevels)
            {
                writer.write(format(",Pct%dx", coverage));
            }

            writer.newLine();

            final Statistics statistics = metrics.statistics();

            int genomeTerritory = metrics.zeroCoverageBases() + metrics.coverageBases();

            writer.write(format("%d,%d,%.3f,%.3f,%.3f",
                    genomeTerritory, metrics.coverageBases(), statistics.Mean, statistics.Median, statistics.StandardDeviation));

            for(FilterType type : FilterType.values())
            {
                if(type == FilterType.UNFILTERED)
                    continue;

                writer.write(format(",%.5f", metrics.calcFilteredPercentage(type)));
            }

            for(Integer coverage : coverageLevels)
            {
                writer.write(format(",%.5f", metrics.calcCoverageFrequency(coverage)));
            }

            writer.newLine();

            writer.close();

            // write coverage frequency for unfiltered aligned bases
            filename = mConfig.formFilename("coverage");
            writer = createBufferedWriter(filename, false);

            writer.write("Coverage,Frequency");
            writer.newLine();

            for(int i = 0; i < metrics.CoverageFrequency.length; ++i)
            {
                writer.write(format("%d,%d", i, metrics.CoverageFrequency[i]));
                writer.newLine();
            }

            writer.close();
        }
        catch(IOException e)
        {
            BM_LOGGER.error("failed to write output file: {}", e.toString());
        }
    }

    public static void main(@NotNull final String[] args) throws Exception
    {
        final VersionInfo version = new VersionInfo("bam-metrics.version");
        BM_LOGGER.info("BamMetrics version: {}", version.version());

        final Options options = createCmdLineOptions();

        try
        {
            final CommandLine cmd = createCommandLine(args, options);

            setLogLevel(cmd);

            BamMetricsApplication bamMetrics = new BamMetricsApplication(cmd);
            bamMetrics.run();
        }
        catch(ParseException e)
        {
            BM_LOGGER.warn(e);
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
