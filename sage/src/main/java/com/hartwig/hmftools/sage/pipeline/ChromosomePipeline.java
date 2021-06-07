package com.hartwig.hmftools.sage.pipeline;

import java.io.File;
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
import com.hartwig.hmftools.common.genome.chromosome.MitochondrialChromosome;
import com.hartwig.hmftools.common.genome.region.GenomeRegion;
import com.hartwig.hmftools.common.variant.hotspot.VariantHotspot;
import com.hartwig.hmftools.sage.config.SageConfig;
import com.hartwig.hmftools.sage.coverage.Coverage;
import com.hartwig.hmftools.sage.phase.Phase;
import com.hartwig.hmftools.sage.quality.QualityRecalibrationMap;
import com.hartwig.hmftools.sage.read.ReadContextCounter;
import com.hartwig.hmftools.sage.variant.SageVariant;
import com.hartwig.hmftools.sage.variant.SageVariantContextFactory;
import com.hartwig.hmftools.sage.variant.SageVariantTier;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import htsjdk.variant.variantcontext.VariantContext;

public class ChromosomePipeline implements AutoCloseable
{

    private static final Logger LOGGER = LogManager.getLogger(ChromosomePipeline.class);
    private static final EnumSet<SageVariantTier> PANEL_ONLY_TIERS = EnumSet.of(SageVariantTier.HOTSPOT, SageVariantTier.PANEL);

    private final String chromosome;
    private final SageConfig config;
    private final List<RegionFuture<List<SageVariant>>> regions = Lists.newArrayList();
    private final IndexedFastaSequenceFile refGenome;
    private final SageVariantPipeline sageVariantPipeline;
    private final Consumer<VariantContext> consumer;
    private final ChromosomePartition partition;
    private final Phase phase;

    public ChromosomePipeline(@NotNull final String chromosome, @NotNull final SageConfig config, @NotNull final Executor executor,
            @NotNull final List<VariantHotspot> hotspots, @NotNull final List<GenomeRegion> panelRegions,
            @NotNull final List<GenomeRegion> highConfidenceRegions, final Map<String, QualityRecalibrationMap> qualityRecalibrationMap,
            @NotNull final Coverage coverage, final Consumer<VariantContext> consumer) throws IOException
    {
        this.chromosome = chromosome;
        this.config = config;
        this.refGenome = new IndexedFastaSequenceFile(new File(config.refGenome()));
        this.consumer = consumer;
        this.sageVariantPipeline = new SomaticPipeline(config,
                executor,
                refGenome,
                hotspots,
                panelRegions,
                highConfidenceRegions,
                qualityRecalibrationMap,
                coverage);
        this.partition = new ChromosomePartition(config, refGenome);
        this.phase = new Phase(config, chromosome, this::write);
    }

    @NotNull
    public String chromosome()
    {
        return chromosome;
    }

    public void process() throws ExecutionException, InterruptedException
    {
        for(GenomeRegion region : partition.partition(chromosome))
        {
            final CompletableFuture<List<SageVariant>> future = sageVariantPipeline.variants(region);
            final RegionFuture<List<SageVariant>> regionFuture = new RegionFuture<>(region, future);
            regions.add(regionFuture);
        }

        submit().get();
    }

    public void process(int minPosition, int maxPosition) throws ExecutionException, InterruptedException
    {
        for(GenomeRegion region : partition.partition(chromosome, minPosition, maxPosition))
        {
            final CompletableFuture<List<SageVariant>> future = sageVariantPipeline.variants(region);
            final RegionFuture<List<SageVariant>> regionFuture = new RegionFuture<>(region, future);
            regions.add(regionFuture);
        }

        submit().get();
    }

    @NotNull
    private CompletableFuture<ChromosomePipeline> submit()
    {
        // Even if regions were executed out of order, they must be phased in order
        regions.sort(Comparator.comparing(RegionFuture::region));

        // Phasing must be done in order but we can do it eagerly as each new region comes in.
        // It is not necessary to wait for the entire chromosome to be finished to start.
        CompletableFuture<Void> done = CompletableFuture.completedFuture(null);
        final Iterator<RegionFuture<List<SageVariant>>> regionsIterator = regions.iterator();
        while(regionsIterator.hasNext())
        {
            CompletableFuture<List<SageVariant>> region = regionsIterator.next().future();
            done = done.thenCombine(region, (aVoid, sageVariants) ->
            {

                sageVariants.forEach(phase);
                return null;
            });

            regionsIterator.remove();
        }

        return done.thenApply(aVoid ->
        {
            phase.flush();
            LOGGER.info("Processing chromosome {} complete", chromosome);
            return ChromosomePipeline.this;
        });
    }

    private void write(@NotNull final SageVariant entry)
    {
        if(include(entry, this.phase.passingPhaseSets()))
        {
            consumer.accept(SageVariantContextFactory.create(entry));
        }
    }

    private boolean include(@NotNull final SageVariant entry, @NotNull final Set<Integer> passingPhaseSets)
    {
        if(config.panelOnly() && !PANEL_ONLY_TIERS.contains(entry.tier()))
        {
            return false;
        }

        if(entry.isPassing())
        {
            return true;
        }

        if(config.filter().hardFilter())
        {
            return false;
        }

        if(entry.tier() == SageVariantTier.HOTSPOT)
        {
            return true;
        }

        // Its not always 100% transparent whats happening with the mixed germline dedup logic unless we keep all the associated records.
        if(entry.mixedGermlineImpact() > 0)
        {
            return true;
        }

        if(!entry.isNormalEmpty() && !entry.isTumorEmpty() && !MitochondrialChromosome.contains(entry.chromosome())
                && !passingPhaseSets.contains(entry.localPhaseSet()))
        {
            final ReadContextCounter normal = entry.normalAltContexts().get(0);
            if(normal.altSupport() > config.filter().filteredMaxNormalAltSupport())
            {
                return false;
            }
        }

        return true;
    }

    @Override
    public void close() throws IOException
    {
        refGenome.close();
    }

    private static class RegionFuture<T>
    {

        private final CompletableFuture<T> future;
        private final GenomeRegion region;

        public RegionFuture(final GenomeRegion region, final CompletableFuture<T> future)
        {
            this.region = region;
            this.future = future;
        }

        public CompletableFuture<T> future()
        {
            return future;
        }

        public GenomeRegion region()
        {
            return region;
        }
    }
}
