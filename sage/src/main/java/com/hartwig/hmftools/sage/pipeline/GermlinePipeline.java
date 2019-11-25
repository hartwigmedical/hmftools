package com.hartwig.hmftools.sage.pipeline;

import java.util.Collections;
import java.util.List;
import java.util.concurrent.CompletableFuture;
import java.util.concurrent.Executor;
import java.util.function.Supplier;
import java.util.stream.Collectors;

import com.hartwig.hmftools.common.genome.region.GenomeRegion;
import com.hartwig.hmftools.common.variant.hotspot.VariantHotspot;
import com.hartwig.hmftools.sage.config.SageConfig;
import com.hartwig.hmftools.sage.context.AltContext;
import com.hartwig.hmftools.sage.context.AltContextSupplier;
import com.hartwig.hmftools.sage.context.RefSequence;
import com.hartwig.hmftools.sage.sam.SamSlicerFactory;
import com.hartwig.hmftools.sage.variant.SageVariant;
import com.hartwig.hmftools.sage.variant.SageVariantFactory;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

import htsjdk.samtools.reference.IndexedFastaSequenceFile;

class GermlinePipeline implements Supplier<CompletableFuture<List<SageVariant>>> {

    private static final Logger LOGGER = LogManager.getLogger(GermlinePipeline.class);

    private final GenomeRegion region;
    private final SageConfig config;
    private final Executor executor;
    private final RefSequence refSequence;
    private final SageVariantFactory variantFactory;
    private final SamSlicerFactory samSlicerFactory;
    private final List<VariantHotspot> hotspots;
    private final List<GenomeRegion> panelRegions;

    GermlinePipeline(final GenomeRegion region, final SageConfig config, final Executor executor, final IndexedFastaSequenceFile refGenome,
            final SamSlicerFactory samSlicerFactory, final List<VariantHotspot> hotspots, final List<GenomeRegion> panelRegions) {
        this.region = region;
        this.config = config;
        this.executor = executor;
        this.variantFactory = new SageVariantFactory(config.filter(), hotspots, panelRegions);
        this.samSlicerFactory = samSlicerFactory;
        this.hotspots = hotspots;
        this.panelRegions = panelRegions;
        this.refSequence = new RefSequence(region, refGenome);
    }

    @NotNull
    public CompletableFuture<List<SageVariant>> submit() {

        final CompletableFuture<List<AltContext>> candidates = CompletableFuture.supplyAsync(new AltContextSupplier(config,
                config.reference(),
                region,
                config.referenceBam(),
                refSequence,
                samSlicerFactory,
                hotspots,
                panelRegions), executor);

        return candidates.thenApply(aVoid -> candidates.join()
                .stream()
                .map(x -> variantFactory.create(x, Collections.emptyList()))
                .collect(Collectors.toList()));
    }

    @Override
    public CompletableFuture<List<SageVariant>> get() {
        return submit();
    }
}
