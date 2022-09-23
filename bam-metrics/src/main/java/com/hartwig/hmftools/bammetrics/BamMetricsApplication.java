package com.hartwig.hmftools.bammetrics;

import static java.lang.String.format;

import static com.hartwig.hmftools.bammetrics.BmConfig.BM_LOGGER;
import static com.hartwig.hmftools.bammetrics.BmConfig.createCmdLineOptions;
import static com.hartwig.hmftools.common.utils.ConfigUtils.setLogLevel;

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

        Metrics combinedMetrics = new Metrics(mConfig.PartitionSize);

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

        long timeTakenMs = System.currentTimeMillis() - startTimeMs;
        double timeTakeMins = timeTakenMs / 60000.0;

        if(mConfig.PerfDebug)
        {
            // combinedStats.PerfCounters.forEach(x -> x.logStats());
        }

        BM_LOGGER.info("BamMetrics complete, mins({})", format("%.3f", timeTakeMins));
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
