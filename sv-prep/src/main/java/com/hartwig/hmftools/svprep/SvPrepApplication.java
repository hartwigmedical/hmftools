package com.hartwig.hmftools.svprep;

import static java.lang.Math.max;
import static java.lang.String.format;

import static com.hartwig.hmftools.common.utils.ConfigUtils.setLogLevel;
import static com.hartwig.hmftools.svprep.SvCommon.SV_LOGGER;
import static com.hartwig.hmftools.svprep.SvConfig.createCmdLineOptions;

import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.utils.version.VersionInfo;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.jetbrains.annotations.NotNull;

public class SvPrepApplication
{
    private final SvConfig mConfig;
    private final ResultsWriter mWriter;
    private final SpanningReadCache mSpanningReadCache;
    private final ExistingJunctionCache mExistingJunctionCache;

    public SvPrepApplication(final CommandLine cmd)
    {
        mConfig = new SvConfig(cmd);
        mWriter = new ResultsWriter(mConfig);
        mSpanningReadCache = new SpanningReadCache(mConfig);
        mExistingJunctionCache = new ExistingJunctionCache();
    }

    public void run()
    {
        if(!mConfig.isValid())
            System.exit(1);

        SV_LOGGER.info("sample({}) starting SV prep", mConfig.SampleId);

        long startTimeMs = System.currentTimeMillis();

        if(mConfig.CalcFragmentLength)
            calcFragmentDistribution();

        mExistingJunctionCache.loadJunctions(mConfig.ExistingJunctionFile);

        CombinedStats combinedStats = new CombinedStats();

        for(HumanChromosome chromosome : HumanChromosome.values())
        {
            String chromosomeStr = mConfig.RefGenVersion.versionedChromosome(chromosome.toString());

            if(!mConfig.SpecificChromosomes.isEmpty() && !mConfig.SpecificChromosomes.contains(chromosomeStr))
                continue;

            SV_LOGGER.info("processing chromosome({})", chromosomeStr);

            ChromosomeTask chromosomeTask = new ChromosomeTask(chromosomeStr, mConfig, mSpanningReadCache, mExistingJunctionCache, mWriter);
            chromosomeTask.process();
            combinedStats.addPartitionStats(chromosomeTask.combinedStats().ReadStats);

            if(combinedStats.PerfCounters.isEmpty())
            {
                combinedStats.PerfCounters.addAll(chromosomeTask.combinedStats().PerfCounters);
            }
            else
            {
                for(int i = 0; i < chromosomeTask.combinedStats().PerfCounters.size(); ++i)
                {
                    combinedStats.PerfCounters.get(i).merge(chromosomeTask.combinedStats().PerfCounters.get((i)));
                }
            }

            System.gc();
        }

        if(mConfig.UseCacheBam)
            mSpanningReadCache.candidateBamWriter().assignCandidateReads(mWriter);

        mWriter.close();

        long timeTakenMs = System.currentTimeMillis() - startTimeMs;
        double timeTakeMins = timeTakenMs / 60000.0;

        if(combinedStats.ReadStats.TotalReads > 10000 || timeTakenMs > 10000)
        {
            SV_LOGGER.info("final stats: {}", combinedStats.ReadStats.toString());
            combinedStats.PerfCounters.forEach(x -> x.logStats());
        }

        SV_LOGGER.info("SvPrep complete, mins({})", format("%.3f", timeTakeMins));
    }

    private void calcFragmentDistribution()
    {
        FragmentSizeDistribution fragSizeDistribution = new FragmentSizeDistribution(mConfig);
        fragSizeDistribution.run();

        final int[] lengthRange = fragSizeDistribution.calcFragmentLengthRange();
        int fragLengthFilterMin = max(mConfig.ReadLength, (int)(lengthRange[0] * 0.5));
        int fragLengthFilterMax = (int)(lengthRange[1] * 1.5);
        mConfig.ReadFiltering.config().setFragmentLengthMin(fragLengthFilterMin, fragLengthFilterMax);
    }

    public static void main(@NotNull final String[] args) throws Exception
    {
        final VersionInfo version = new VersionInfo("sv-prep.version");
        SV_LOGGER.info("sv-prep version: {}", version.version());

        final Options options = createCmdLineOptions();

        try
        {
            final CommandLine cmd = createCommandLine(args, options);

            setLogLevel(cmd);

            SvPrepApplication svPrep = new SvPrepApplication(cmd);
            svPrep.run();
        }
        catch(ParseException e)
        {
            SV_LOGGER.warn(e);
            final HelpFormatter formatter = new HelpFormatter();
            formatter.printHelp("Isofox", options);
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
