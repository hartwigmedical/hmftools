package com.hartwig.hmftools.sage.pipeline;

import java.util.List;
import java.util.concurrent.CompletableFuture;
import java.util.function.Consumer;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.sage.config.SageConfig;
import com.hartwig.hmftools.sage.context.AltContext;
import com.hartwig.hmftools.sage.phase.Phase;
import com.hartwig.hmftools.sage.variant.SageVariant;
import com.hartwig.hmftools.sage.variant.SageVariantContextFactory;
import com.hartwig.hmftools.sage.variant.SageVariantFactory;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import htsjdk.variant.variantcontext.VariantContext;

public class ChromosomePipeline implements Consumer<CompletableFuture<List<SageVariant>>> {

    private static final Logger LOGGER = LogManager.getLogger(ChromosomePipeline.class);

    private final String chromosome;
    private final SageConfig config;
    private final IndexedFastaSequenceFile reference;
    private final SageVariantFactory sageVariantFactory;
    private final List<VariantContext> variantContexts = Lists.newArrayList();
    private final List<CompletableFuture<List<SageVariant>>> regions = Lists.newArrayList();

    public ChromosomePipeline(@NotNull final String chromosome, @NotNull final SageConfig config,
            @NotNull final IndexedFastaSequenceFile reference, @NotNull final SageVariantFactory sageVariantFactory) {
        this.chromosome = chromosome;
        this.config = config;
        this.reference = reference;
        this.sageVariantFactory = sageVariantFactory;
    }

    @NotNull
    public String chromosome() {
        return chromosome;
    }

    @Override
    public void accept(final CompletableFuture<List<SageVariant>> regionData) {
        regions.add(regionData);
    }

    @NotNull
    public List<VariantContext> variantContexts() {
        return variantContexts;
    }

    @NotNull
    public CompletableFuture<ChromosomePipeline> submit() {

        final CompletableFuture<Void> doneChromosome = CompletableFuture.allOf(regions.toArray(new CompletableFuture[regions.size()]));
        return doneChromosome.thenApply((aVoid) -> {

            LOGGER.info("Phasing chromosome {}", chromosome);

            final Consumer<SageVariant> phasedConsumer = variant -> {
                if (include(variant)) {
                    final VariantContext context = SageVariantContextFactory.create(variant);
                    variantContexts.add(context);
                }
            };

            final Phase phase = new Phase(reference, sageVariantFactory, phasedConsumer);
            for (CompletableFuture<List<SageVariant>> region : regions) {
                for (SageVariant sageVariant : region.join()) {
                    phase.accept(sageVariant);
                }
            }

            phase.flush();
            return this;
        });
    }

    private boolean include(@NotNull final SageVariant entry) {
        if (entry.isPassing()) {
            return true;
        }

        if (config.filter().hardFilter()) {
            return false;
        }

        final AltContext normal = entry.normal();
        return normal.altSupport() <= config.filter().hardMaxNormalAltSupport();
    }

}
