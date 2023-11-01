package com.hartwig.hmftools.svprep;

import static java.lang.String.format;

import static com.hartwig.hmftools.common.utils.PerformanceCounter.runTimeMinsStr;
import static com.hartwig.hmftools.common.utils.config.ConfigUtils.setLogLevel;
import static com.hartwig.hmftools.svprep.SvCommon.APP_NAME;
import static com.hartwig.hmftools.svprep.SvCommon.SV_LOGGER;

import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;
import com.hartwig.hmftools.common.utils.version.VersionInfo;

import org.jetbrains.annotations.NotNull;

public class SvPrepApplication
{
    private final SvConfig mConfig;
    private final ResultsWriter mWriter;
    private final SpanningReadCache mSpanningReadCache;
    private final ExistingJunctionCache mExistingJunctionCache;

    public SvPrepApplication(final ConfigBuilder configBuilder)
    {
        mConfig = new SvConfig(configBuilder);
        mWriter = new ResultsWriter(mConfig);
        mSpanningReadCache = new SpanningReadCache(mConfig);
        mExistingJunctionCache = new ExistingJunctionCache();
    }

    public void run()
    {
        if(!mConfig.isValid())
            System.exit(1);

        SV_LOGGER.info("running SvPrep for sample({})", mConfig.SampleId);

        long startTimeMs = System.currentTimeMillis();

        if(mConfig.CalcFragmentLength)
            calcFragmentDistribution();

        mExistingJunctionCache.loadJunctions(mConfig.ExistingJunctionFile);

        CombinedStats combinedStats = new CombinedStats();

        for(HumanChromosome chromosome : HumanChromosome.values())
        {
            String chromosomeStr = mConfig.RefGenVersion.versionedChromosome(chromosome.toString());

            if(mConfig.SpecificChrRegions.excludeChromosome(chromosomeStr))
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
        }

        if(mConfig.UseCacheBam)
        {
            mSpanningReadCache.reset(); // clear data before candidate assignment
            mSpanningReadCache.candidateBamWriter().assignCandidateReads(mWriter);
        }

        mWriter.close();

        long timeTakenMs = System.currentTimeMillis() - startTimeMs;

        if(mConfig.PerfDebug && (combinedStats.ReadStats.TotalReads > 10000 || timeTakenMs > 10000))
        {
            SV_LOGGER.info("final stats: {}", combinedStats.ReadStats.toString());

            if(SV_LOGGER.isDebugEnabled())
                combinedStats.PerfCounters.forEach(x -> x.logIntervalStats());
            else
                combinedStats.PerfCounters.forEach(x -> x.logStats());
        }

        SV_LOGGER.info("SvPrep complete, mins({})", runTimeMinsStr(startTimeMs));
    }

    private void calcFragmentDistribution()
    {
        FragmentSizeDistribution fragSizeDistribution = new FragmentSizeDistribution(mConfig);
        fragSizeDistribution.run();

        final int[] fragmentDistLengths = fragSizeDistribution.calculatePercentileLengths();

        SV_LOGGER.info("fragment length distribution percentile values(min={} max={})",
                fragmentDistLengths[0], fragmentDistLengths[1]);

        if(fragmentDistLengths[0] > 0 && fragmentDistLengths[1] > 0)
            mConfig.ReadFiltering.config().setFragmentLengths(fragmentDistLengths[0], fragmentDistLengths[1]);

        if(fragSizeDistribution.maxReadLength() > 0)
            mConfig.ReadLength = fragSizeDistribution.maxReadLength();
    }

    public static void main(@NotNull final String[] args)
    {
        ConfigBuilder configBuilder = new ConfigBuilder(APP_NAME);
        SvConfig.addConfig(configBuilder);

        configBuilder.checkAndParseCommandLine(args);

        SvPrepApplication svPrep = new SvPrepApplication(configBuilder);
        svPrep.run();
    }
}
