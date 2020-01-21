package com.hartwig.hmftools.sage.pipeline;

import java.util.List;
import java.util.concurrent.CompletableFuture;
import java.util.concurrent.Executor;
import java.util.stream.Collectors;

import com.hartwig.hmftools.common.genome.region.GenomeRegion;
import com.hartwig.hmftools.common.variant.hotspot.VariantHotspot;
import com.hartwig.hmftools.sage.config.SageConfig;
import com.hartwig.hmftools.sage.context.AltContext;
import com.hartwig.hmftools.sage.context.RefSequence;
import com.hartwig.hmftools.sage.evidence.PrimaryEvidence;
import com.hartwig.hmftools.sage.sam.SamSlicerFactory;
import com.hartwig.hmftools.sage.variant.SageVariant;
import com.hartwig.hmftools.sage.variant.SageVariantFactory;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

class GermlineOnlyPipeline implements SageVariantPipeline {

    private static final Logger LOGGER = LogManager.getLogger(GermlineOnlyPipeline.class);

    private final SageConfig config;
    private final Executor executor;
    private final List<VariantHotspot> hotspots;
    private final List<GenomeRegion> panelRegions;
    private final PrimaryEvidence primaryEvidence;

    GermlineOnlyPipeline(final SageConfig config, final Executor executor, final List<VariantHotspot> hotspots,
            final List<GenomeRegion> panelRegions) {
        this.config = config;
        this.executor = executor;
        this.hotspots = hotspots;
        this.panelRegions = panelRegions;

        final SamSlicerFactory samSlicerFactory = new SamSlicerFactory(config, panelRegions);
        this.primaryEvidence = new PrimaryEvidence(config, hotspots, panelRegions, samSlicerFactory);

    }

    @NotNull
    @Override
    public CompletableFuture<List<SageVariant>> variants(@NotNull final GenomeRegion region, @NotNull final RefSequence refSequence) {

        final SageVariantFactory variantFactory = new SageVariantFactory(config.filter(), hotspots, panelRegions);
        final CompletableFuture<List<AltContext>> candidates =
                CompletableFuture.supplyAsync(() -> primaryEvidence.get(config.reference(), config.referenceBam(), refSequence, region),
                        executor);

        return candidates.thenApply(aVoid -> candidates.join().stream().map(variantFactory::create).collect(Collectors.toList()));
    }

}
