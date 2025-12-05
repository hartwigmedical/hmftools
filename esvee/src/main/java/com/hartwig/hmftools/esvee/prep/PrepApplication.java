package com.hartwig.hmftools.esvee.prep;

import static com.hartwig.hmftools.common.perf.PerformanceCounter.runTimeMinsStr;
import static com.hartwig.hmftools.esvee.assembly.AssemblyConfig.SV_LOGGER;
import static com.hartwig.hmftools.esvee.common.FileCommon.APP_NAME;
import static com.hartwig.hmftools.esvee.prep.types.DiscordantStats.writeDiscordantStats;

import java.util.stream.Collectors;

import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;
import com.hartwig.hmftools.esvee.common.FragmentLengthBounds;
import com.hartwig.hmftools.esvee.prep.types.CombinedStats;

import org.jetbrains.annotations.NotNull;

public class PrepApplication
{
    private final PrepConfig mConfig;
    private final ResultsWriter mWriter;
    private final SpanningReadCache mSpanningReadCache;

    public PrepApplication(final ConfigBuilder configBuilder)
    {
        mConfig = new PrepConfig(configBuilder);
        mWriter = new ResultsWriter(mConfig);
        mSpanningReadCache = new SpanningReadCache(mConfig);
    }

    public void run()
    {
        if(!mConfig.isValid())
            System.exit(1);

        SV_LOGGER.info("running Esvee Prep for {}",
                mConfig.SampleIds.size() == 1 ? mConfig.sampleId() : mConfig.SampleIds.stream().collect(Collectors.joining(",")));

        long startTimeMs = System.currentTimeMillis();

        calcFragmentDistribution();

        CombinedStats combinedStats = new CombinedStats();

        for(HumanChromosome chromosome : HumanChromosome.values())
        {
            String chromosomeStr = mConfig.RefGenVersion.versionedChromosome(chromosome.toString());

            if(mConfig.SpecificChrRegions.excludeChromosome(chromosomeStr))
                continue;

            SV_LOGGER.info("processing chromosome({})", chromosomeStr);

            ChromosomeTask chromosomeTask = new ChromosomeTask(chromosomeStr, mConfig, mSpanningReadCache, mWriter);
            chromosomeTask.process();
            combinedStats.addPartitionStats(chromosomeTask.combinedStats().ReadStats);
            combinedStats.addDiscordantStats(chromosomeTask.combinedStats().Discordants);

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

        writeDiscordantStats(mConfig, combinedStats.ReadStats.TotalReads, mWriter.writtenCount(), combinedStats.Discordants);

        if(mConfig.UseCacheBam)
        {
            mSpanningReadCache.reset(); // clear data before candidate assignment
            mSpanningReadCache.candidateBamWriter().assignCandidateReads(mWriter);
        }

        if(mConfig.BamToolPath != null)
            System.gc(); // call to free up prep memory before the BAM tools sort & index routines are run

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

        SV_LOGGER.info("Esvee prep complete, mins({})", runTimeMinsStr(startTimeMs));
    }

    private void calcFragmentDistribution()
    {
        FragmentSizeDistribution fragSizeDistribution = new FragmentSizeDistribution(mConfig);
        fragSizeDistribution.run();

        FragmentLengthBounds fragmentLengthBounds = fragSizeDistribution.calculateFragmentLengthBounds();

        SV_LOGGER.info("fragment length distribution bounds(min={} max={} median={})",
                fragmentLengthBounds.LowerBound, fragmentLengthBounds.UpperBound, fragmentLengthBounds.Median);

        if(fragmentLengthBounds.isValid())
        {
            mConfig.ReadFiltering.config().setFragmentLengths(
                    fragmentLengthBounds.LowerBound,
                    mConfig.MaxFragmentLengthOverride > 0 ? mConfig.MaxFragmentLengthOverride : fragmentLengthBounds.UpperBound);
        }
    }

    public static void main(@NotNull final String[] args)
    {
        ConfigBuilder configBuilder = new ConfigBuilder(APP_NAME);
        PrepConfig.registerConfig(configBuilder);

        configBuilder.checkAndParseCommandLine(args);

        PrepApplication svPrep = new PrepApplication(configBuilder);
        svPrep.run();
    }
}
