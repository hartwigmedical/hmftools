package com.hartwig.hmftools.sage.pipeline;

import static com.hartwig.hmftools.sage.ReferenceData.loadRefGenome;
import static com.hartwig.hmftools.sage.SageCommon.SG_LOGGER;

import java.io.IOException;
import java.util.EnumSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.concurrent.Callable;
import java.util.concurrent.CompletableFuture;
import java.util.concurrent.Executor;
import java.util.function.Consumer;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.genome.chromosome.Chromosome;
import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.genome.chromosome.MitochondrialChromosome;
import com.hartwig.hmftools.common.utils.PerformanceCounter;
import com.hartwig.hmftools.common.utils.TaskExecutor;
import com.hartwig.hmftools.common.utils.sv.ChrBaseRegion;
import com.hartwig.hmftools.sage.ReferenceData;
import com.hartwig.hmftools.sage.config.SageConfig;
import com.hartwig.hmftools.sage.coverage.Coverage;
import com.hartwig.hmftools.sage.phase.PhaseSetCounter;
import com.hartwig.hmftools.sage.phase.VariantDeduper;
import com.hartwig.hmftools.sage.quality.QualityRecalibrationMap;
import com.hartwig.hmftools.sage.evidence.ReadContextCounter;
import com.hartwig.hmftools.sage.common.SageVariant;
import com.hartwig.hmftools.sage.common.VariantTier;

import htsjdk.samtools.reference.IndexedFastaSequenceFile;

public class ChromosomePipeline implements AutoCloseable
{
    private final String mChromosome;
    private final SageConfig mConfig;
    private final List<RegionFuture<List<SageVariant>>> mRegions = Lists.newArrayList();
    private final IndexedFastaSequenceFile mRefGenome;
    private final Consumer<SageVariant> mWriteConsumer;
    private final ChromosomePartition mPartition;
    private final List<RegionTask> mRegionTasks;

    private static final EnumSet<VariantTier> PANEL_ONLY_TIERS = EnumSet.of(VariantTier.HOTSPOT, VariantTier.PANEL);

    public ChromosomePipeline(
            final String chromosome, final SageConfig config,
            final ReferenceData refData, final Map<String,QualityRecalibrationMap> qualityRecalibrationMap,
            final Coverage coverage, final PhaseSetCounter phaseSetCounter, final Consumer<SageVariant> consumer)
    {
        mChromosome = chromosome;
        mConfig = config;
        mRefGenome = loadRefGenome(config.RefGenomeFile);
        mWriteConsumer = consumer;

        final Chromosome chr = HumanChromosome.contains(chromosome)
                ? HumanChromosome.fromString(chromosome) : MitochondrialChromosome.fromString(chromosome);

        mPartition = new ChromosomePartition(config, mRefGenome);

        List<ChrBaseRegion> partitionedRegions = mPartition.partition(mChromosome);

        mRegionTasks = Lists.newArrayList();

        for(int i = 0; i < partitionedRegions.size(); ++i)
        {
            ChrBaseRegion region = partitionedRegions.get(i);

            mRegionTasks.add(new RegionTask(i, region, config, mRefGenome,
                    refData.Hotspots.get(chr), refData.PanelWithHotspots.get(chr), refData.ChromosomeTranscripts.get(chromosome),
                    refData.HighConfidence.get(chr), qualityRecalibrationMap, phaseSetCounter, coverage));
        }
    }

    public String chromosome()
    {
        return mChromosome;
    }

    public void process()
    {
        SG_LOGGER.info("chromosome({}) executing {} regions", mChromosome, mRegionTasks.size());

        final List<Callable> callableList = mRegionTasks.stream().collect(Collectors.toList());
        TaskExecutor.executeTasks(callableList, mConfig.Threads);

        SG_LOGGER.debug("chromosome({}) {} regions complete", mChromosome, mRegionTasks.size());

        // write variants to file
        mRegionTasks.forEach(x -> x.writeVariants(mWriteConsumer));

        if(SG_LOGGER.isDebugEnabled())
        {
            List<PerformanceCounter> perfCounters = mRegionTasks.get(0).getPerfCounters();

            for(int i = 1; i < mRegionTasks.size(); ++i)
            {
                List<PerformanceCounter> taskPerfCounters = mRegionTasks.get(i).getPerfCounters();

                for(int j = 0; j < perfCounters.size(); ++j)
                {
                    perfCounters.get(j).merge(taskPerfCounters.get(j));
                }
            }

            perfCounters.forEach(x -> x.logStats());
        }

        SG_LOGGER.info("chromosome({}) analysis complete", mChromosome);
    }

    @Override
    public void close() throws IOException
    {
        mRefGenome.close();
    }

    private static class RegionFuture<T>
    {
        private final CompletableFuture<T> mFuture;
        private final ChrBaseRegion mRegion;

        public RegionFuture(final ChrBaseRegion region, final CompletableFuture<T> future)
        {
            mRegion = region;
            mFuture = future;
        }

        public CompletableFuture<T> future()
        {
            return mFuture;
        }

        public ChrBaseRegion region()
        {
            return mRegion;
        }
    }
}
