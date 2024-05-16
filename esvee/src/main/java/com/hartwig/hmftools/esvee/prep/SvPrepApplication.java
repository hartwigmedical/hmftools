package com.hartwig.hmftools.esvee.prep;

import static com.hartwig.hmftools.common.utils.PerformanceCounter.runTimeMinsStr;
import static com.hartwig.hmftools.esvee.AssemblyConfig.SV_LOGGER;
import static com.hartwig.hmftools.esvee.common.FileCommon.APP_NAME;

import java.util.stream.Collectors;

import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;
import com.hartwig.hmftools.esvee.common.FragmentLengthBounds;
import com.hartwig.hmftools.esvee.prep.types.CombinedStats;

import org.jetbrains.annotations.NotNull;

public class SvPrepApplication
{
    private final PrepConfig mConfig;
    private final ResultsWriter mWriter;
    private final SpanningReadCache mSpanningReadCache;

    public SvPrepApplication(final ConfigBuilder configBuilder)
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

        if(mConfig.CalcFragmentLength)
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

        SV_LOGGER.info("Esvee Prep complete, mins({})", runTimeMinsStr(startTimeMs));
    }

    private void calcFragmentDistribution()
    {
        FragmentSizeDistribution fragSizeDistribution = new FragmentSizeDistribution(mConfig);
        fragSizeDistribution.run();

        FragmentLengthBounds fragmentLengthBounds = fragSizeDistribution.calculateFragmentLengthBounds();

        SV_LOGGER.info("fragment length distribution bounds(min={} max={})",
                fragmentLengthBounds.LowerBound, fragmentLengthBounds.UpperBound);

        if(fragmentLengthBounds.isValid())
            mConfig.ReadFiltering.config().setFragmentLengths(fragmentLengthBounds.LowerBound, fragmentLengthBounds.UpperBound);

        if(fragSizeDistribution.maxReadLength() > 0)
            mConfig.ReadLength = fragSizeDistribution.maxReadLength();
    }

    public static void main(@NotNull final String[] args)
    {
        ConfigBuilder configBuilder = new ConfigBuilder(APP_NAME);
        PrepConfig.addConfig(configBuilder);

        configBuilder.checkAndParseCommandLine(args);

        SvPrepApplication svPrep = new SvPrepApplication(configBuilder);
        svPrep.run();
    }
}
