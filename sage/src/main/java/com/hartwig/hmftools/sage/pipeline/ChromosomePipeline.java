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
import java.util.concurrent.CompletableFuture;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.Executor;
import java.util.function.Consumer;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.genome.chromosome.Chromosome;
import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.genome.chromosome.MitochondrialChromosome;
import com.hartwig.hmftools.common.utils.sv.ChrBaseRegion;
import com.hartwig.hmftools.sage.ReferenceData;
import com.hartwig.hmftools.sage.config.SageConfig;
import com.hartwig.hmftools.sage.coverage.Coverage;
import com.hartwig.hmftools.sage.phase.PhaseSetCounter;
import com.hartwig.hmftools.sage.phase.VariantPhaser;
import com.hartwig.hmftools.sage.quality.QualityRecalibrationMap;
import com.hartwig.hmftools.sage.evidence.ReadContextCounter;
import com.hartwig.hmftools.sage.variant.SageVariant;
import com.hartwig.hmftools.sage.variant.VariantTier;

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
    private final VariantPhaser mVariantPhaser;

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
                refData.HighConfidence.get(chr), qualityRecalibrationMap, coverage);

        mPartition = new ChromosomePartition(config, mRefGenome);

        mVariantPhaser = new VariantPhaser(refData.ChromosomeTranscripts.get(chromosome), phaseSetCounter, this::write);
    }

    public String chromosome()
    {
        return mChromosome;
    }

    public void process() throws ExecutionException, InterruptedException
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
            ChrBaseRegion region = regionsIterator.next().region();
            CompletableFuture<List<SageVariant>> regionVariantsFuture = regionsIterator.next().future();

            done = done.thenCombine(regionVariantsFuture, (aVoid, sageVariants) ->
            {
                SG_LOGGER.trace("region({}) phasing {} variants", region, sageVariants.size());
                sageVariants.forEach(mVariantPhaser);
                return null;
            });

            regionsIterator.remove();
        }

        return done.thenApply(aVoid ->
        {
            mVariantPhaser.flush();
            SG_LOGGER.info("processing chromosome {} complete", mChromosome);
            return ChromosomePipeline.this;
        });
    }

    private void write(final SageVariant variant)
    {
        if(checkWriteVariant(variant, mVariantPhaser.passingPhaseSets()))
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
