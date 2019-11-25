package com.hartwig.hmftools.sage.pipeline;

import java.io.IOException;
import java.util.List;
import java.util.concurrent.CompletableFuture;
import java.util.concurrent.Executor;
import java.util.function.Function;

import com.google.common.collect.ListMultimap;
import com.hartwig.hmftools.common.genome.chromosome.Chromosome;
import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.genome.region.GenomeRegion;
import com.hartwig.hmftools.common.variant.hotspot.VariantHotspot;
import com.hartwig.hmftools.sage.config.SageConfig;
import com.hartwig.hmftools.sage.sam.SamSlicerFactory;
import com.hartwig.hmftools.sage.variant.SageVariant;
import com.hartwig.hmftools.sage.variant.SageVariantContextFactory;
import com.hartwig.hmftools.sage.variant.SageVariantFactory;

import org.jetbrains.annotations.NotNull;

import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import htsjdk.variant.variantcontext.VariantContext;

public class PipelineFactory {
    private final SageConfig config;
    private final Executor executor;
    private final IndexedFastaSequenceFile refSequence;
    private final SamSlicerFactory samSlicerFactory;
    private final Function<SageVariant, VariantContext> variantContextFactory;

    private final ListMultimap<Chromosome, GenomeRegion> panel;
    private final ListMultimap<Chromosome, VariantHotspot> hotspots;

    public PipelineFactory(@NotNull final SageConfig config, @NotNull final Executor executor,
            @NotNull final IndexedFastaSequenceFile refSequence, @NotNull final SamSlicerFactory samSlicerFactory,
            @NotNull final ListMultimap<Chromosome, VariantHotspot> hotspots, @NotNull final ListMultimap<Chromosome, GenomeRegion> panel) {
        this.config = config;
        this.executor = executor;
        this.refSequence = refSequence;
        this.samSlicerFactory = samSlicerFactory;
        this.hotspots = hotspots;
        this.panel = panel;
        variantContextFactory =
                config.germlineOnly() ? SageVariantContextFactory::germlineOnly : SageVariantContextFactory::pairedTumorNormal;
    }

    public ChromosomePipeline create(@NotNull final String contig) throws IOException {
        final Chromosome chromosome = HumanChromosome.fromString(contig);
        final SageVariantFactory variantFactory = new SageVariantFactory(config.filter(), hotspots.get(chromosome), panel.get(chromosome));

        return new ChromosomePipeline(contig, config, refSequence, variantFactory, variantContextFactory);

    }

    @NotNull
    public CompletableFuture<List<SageVariant>> createPipeline(@NotNull final GenomeRegion region) {
        final Chromosome chromosome = HumanChromosome.fromString(region.chromosome());

        final List<VariantHotspot> chromosomeHotspots = hotspots.get(chromosome);
        final List<GenomeRegion> chromosomePanel = panel.get(chromosome);

        if (config.germlineOnly()) {
            return new GermlinePipeline(region,
                    config,
                    executor,
                    refSequence,
                    samSlicerFactory,
                    chromosomeHotspots,
                    chromosomePanel).submit();
        } else {
            return new SomaticPipeline(region,
                    config,
                    executor,
                    refSequence,
                    samSlicerFactory,
                    chromosomeHotspots,
                    chromosomePanel).submit();
        }
    }
}
