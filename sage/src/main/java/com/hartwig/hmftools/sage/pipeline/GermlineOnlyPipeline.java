package com.hartwig.hmftools.sage.pipeline;

import java.util.List;
import java.util.concurrent.CompletableFuture;
import java.util.concurrent.Executor;
import java.util.function.Predicate;
import java.util.stream.Collectors;

import com.hartwig.hmftools.common.genome.region.GenomeRegion;
import com.hartwig.hmftools.common.variant.hotspot.VariantHotspot;
import com.hartwig.hmftools.sage.config.SageConfig;
import com.hartwig.hmftools.sage.context.AltContext;
import com.hartwig.hmftools.sage.evidence.GermlineEvidence;
import com.hartwig.hmftools.sage.read.ReadContextCounter;
import com.hartwig.hmftools.sage.ref.RefSequence;
import com.hartwig.hmftools.sage.sam.SamSlicerFactory;
import com.hartwig.hmftools.sage.select.HotspotSelector;
import com.hartwig.hmftools.sage.variant.SageVariant;
import com.hartwig.hmftools.sage.variant.SageVariantFactory;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

import htsjdk.samtools.reference.ReferenceSequenceFile;

class GermlineOnlyPipeline implements SageVariantPipeline {

    private static final Logger LOGGER = LogManager.getLogger(GermlineOnlyPipeline.class);

    private final SageConfig config;
    private final Executor executor;
    private final List<VariantHotspot> hotspots;
    private final List<GenomeRegion> panelRegions;
    private final List<GenomeRegion> highConfidenceRegions;
    private final GermlineEvidence germlineEvidence;
    private final ReferenceSequenceFile refGenome;

    GermlineOnlyPipeline(final SageConfig config, final Executor executor, final ReferenceSequenceFile refGenome,
            final List<VariantHotspot> hotspots, final List<GenomeRegion> panelRegions, final List<GenomeRegion> highConfidenceRegions) {
        this.config = config;
        this.executor = executor;
        this.hotspots = hotspots;
        this.panelRegions = panelRegions;

        final SamSlicerFactory samSlicerFactory = new SamSlicerFactory(config, panelRegions);
        this.germlineEvidence = new GermlineEvidence(config, samSlicerFactory, refGenome);
        this.highConfidenceRegions = highConfidenceRegions;
        this.refGenome = refGenome;

    }

    @NotNull
    @Override
    public CompletableFuture<List<SageVariant>> variants(@NotNull final GenomeRegion region) {
        final Predicate<ReadContextCounter> hardFilter = hardFilterEvidence(hotspots);
        final SageVariantFactory variantFactory = new SageVariantFactory(config.filter(), hotspots, panelRegions, highConfidenceRegions);

        final CompletableFuture<RefSequence> refSequenceFuture =
                CompletableFuture.supplyAsync(() -> new RefSequence(region, refGenome), executor);

        final CompletableFuture<List<AltContext>> candidates =
                refSequenceFuture.thenApply(refSequence -> germlineEvidence.get(config.reference().get(0),
                        config.referenceBam().get(0),
                        refSequence,
                        region));

        return candidates.thenApply(aVoid -> candidates.join()
                .stream()
                .map(AltContext::primaryReadContext)
                .filter(hardFilter)
                .map(variantFactory::create)
                .collect(Collectors.toList()));
    }

    @NotNull
    private Predicate<ReadContextCounter> hardFilterEvidence(@NotNull final List<VariantHotspot> variants) {
        return config.filter().hardFilter(new HotspotSelector(variants));
    }

}
