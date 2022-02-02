package com.hartwig.hmftools.sage.pipeline;

import static com.hartwig.hmftools.sage.ReferenceData.loadRefGenome;
import static com.hartwig.hmftools.sage.SageCommon.SG_LOGGER;

import java.io.IOException;
import java.util.Comparator;
import java.util.EnumSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.concurrent.Callable;
import java.util.concurrent.CompletableFuture;
import java.util.concurrent.ExecutionException;
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
    private final SomaticPipeline mSomaticPipeline;
    private final Consumer<SageVariant> mWriteConsumer;
    private final ChromosomePartition mPartition;
    private final List<RegionTask> mRegionTasks;
    private final VariantDeduper mVariantDeduper;

    private static final EnumSet<VariantTier> PANEL_ONLY_TIERS = EnumSet.of(VariantTier.HOTSPOT, VariantTier.PANEL);

    public ChromosomePipeline(
            final String chromosome, final SageConfig config, final Executor executor,
            final ReferenceData refData, final Map<String,QualityRecalibrationMap> qualityRecalibrationMap,
            final Coverage coverage, final PhaseSetCounter phaseSetCounter, final Consumer<SageVariant> consumer)
    {
        mChromosome = chromosome;
        mConfig = config;
        mRefGenome = loadRefGenome(config.RefGenomeFile);
        mWriteConsumer = consumer;

        final Chromosome chr = HumanChromosome.contains(chromosome)
                ? HumanChromosome.fromString(chromosome) : MitochondrialChromosome.fromString(chromosome);

        mSomaticPipeline = new SomaticPipeline(
                config, executor, mRefGenome,
                refData.Hotspots.get(chr), refData.PanelWithHotspots.get(chr),
                refData.HighConfidence.get(chr), qualityRecalibrationMap, phaseSetCounter, coverage);

        mPartition = new ChromosomePartition(config, mRefGenome);

        List<ChrBaseRegion> partitionedRegions = mPartition.partition(mChromosome);

        mRegionTasks = Lists.newArrayList();

        for(int i = 0; i < partitionedRegions.size(); ++i)
        {
            ChrBaseRegion region = partitionedRegions.get(i);

            mRegionTasks.add(new RegionTask(i, region, config, mRefGenome,
                    refData.Hotspots.get(chr), refData.PanelWithHotspots.get(chr),
                    refData.HighConfidence.get(chr), qualityRecalibrationMap, phaseSetCounter, coverage));
        }

        mVariantDeduper = new VariantDeduper(refData.ChromosomeTranscripts.get(chromosome), phaseSetCounter, this::write);
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

        PerformanceCounter perfCounter = new PerformanceCounter("Dedup");

        perfCounter.start();

        for(RegionTask regionTask : mRegionTasks)
        {
            List<SageVariant> regionVariants = regionTask.getVariants();

            SG_LOGGER.trace("phasing {} variants", regionVariants.size());
            regionVariants.forEach(mVariantDeduper);
        }

        mVariantDeduper.flush();

        perfCounter.stop();
        perfCounter.logStats();

        SG_LOGGER.info("chromosome({}) analysis complete", mChromosome);
    }

    public void processOld() throws ExecutionException, InterruptedException
    {
        for(ChrBaseRegion region : mPartition.partition(mChromosome))
        {
            final CompletableFuture<List<SageVariant>> future = mSomaticPipeline.findVariants(region);
            final RegionFuture<List<SageVariant>> regionFuture = new RegionFuture<>(region, future);
            mRegions.add(regionFuture);
        }

        submit().get();
    }

    private CompletableFuture<ChromosomePipeline> submit()
    {
        // even if regions were executed out of order, they must be phased in order
        mRegions.sort(Comparator.comparing(RegionFuture::region));

        // SG_LOGGER.trace("sorted {} completed(?) regions", mRegions.size());

        // Phasing must be done in order but we can do it eagerly as each new region comes in
        // It is not necessary to wait for the entire chromosome to be finished to start
        CompletableFuture<Void> done = CompletableFuture.completedFuture(null);
        final Iterator<RegionFuture<List<SageVariant>>> regionsIterator = mRegions.iterator();
        while(regionsIterator.hasNext())
        {
            // ChrBaseRegion region = regionsIterator.next().region();
            CompletableFuture<List<SageVariant>> regionVariantsFuture = regionsIterator.next().future();

            done = done.thenCombine(regionVariantsFuture, (aVoid, sageVariants) ->
            {
                SG_LOGGER.trace("phasing {} variants", sageVariants.size());
                sageVariants.forEach(mVariantDeduper);
                return null;
            });

            regionsIterator.remove();
        }

        return done.thenApply(aVoid ->
        {
            mVariantDeduper.flush();
            SG_LOGGER.info("processing chromosome {} complete", mChromosome);
            return ChromosomePipeline.this;
        });
    }

    private void write(final SageVariant variant)
    {
        if(checkWriteVariant(variant, mVariantDeduper.passingPhaseSets()))
        {
            mWriteConsumer.accept(variant);
        }
    }

    private boolean checkWriteVariant(final SageVariant variant, final Set<Integer> passingPhaseSets)
    {
        if(mConfig.PanelOnly && !PANEL_ONLY_TIERS.contains(variant.tier()))
            return false;

        if(variant.isPassing())
            return true;

        if(mConfig.Filter.HardFilter)
            return false;

        if(variant.tier() == VariantTier.HOTSPOT)
            return true;

        // Its not always 100% transparent whats happening with the mixed germline dedup logic unless we keep all the associated records
        if(variant.mixedGermlineImpact() > 0)
            return true;

        if(!variant.isNormalEmpty() && !variant.isTumorEmpty() && !MitochondrialChromosome.contains(variant.chromosome())
        && !passingPhaseSets.contains(variant.localPhaseSet()))
        {
            final ReadContextCounter normal = variant.normalAltContexts().get(0);

            if(normal.altSupport() > mConfig.Filter.FilteredMaxNormalAltSupport)
                return false;
        }

        return true;
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
