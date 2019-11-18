package com.hartwig.hmftools.sage.pipeline;

import java.util.List;
import java.util.concurrent.CompletableFuture;
import java.util.function.Consumer;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.sage.SagePipelineData;
import com.hartwig.hmftools.sage.SageVariantContextFactory;
import com.hartwig.hmftools.sage.config.SageConfig;
import com.hartwig.hmftools.sage.context.AltContext;
import com.hartwig.hmftools.sage.phase.Phase;
import com.hartwig.hmftools.sage.variant.SageVariant;
import com.hartwig.hmftools.sage.variant.SageVariantFactory;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import htsjdk.variant.variantcontext.VariantContext;

public class ChromosomePipeline implements Consumer<CompletableFuture<SagePipelineData>> {

    private static final Logger LOGGER = LogManager.getLogger(ChromosomePipeline.class);

    private final String chromosome;
    private final SageConfig config;
    private final IndexedFastaSequenceFile reference;
    private final SageVariantFactory sageVariantFactory;
    private final List<CompletableFuture<SagePipelineData>> regions = Lists.newArrayList();

    public ChromosomePipeline(@NotNull final String chromosome, @NotNull final SageConfig config,
            @NotNull final IndexedFastaSequenceFile reference, @NotNull final SageVariantFactory sageVariantFactory) {
        this.chromosome = chromosome;
        this.config = config;
        this.reference = reference;
        this.sageVariantFactory = sageVariantFactory;
    }

    @Override
    public void accept(final CompletableFuture<SagePipelineData> regionData) {
        regions.add(regionData);
    }

    @NotNull
    public CompletableFuture<List<VariantContext>> submit() {

        final CompletableFuture<Void> doneChromosome = CompletableFuture.allOf(regions.toArray(new CompletableFuture[regions.size()]));
        return doneChromosome.thenApply((Void aVoid) -> {

            LOGGER.info("Phasing chromosome {}", chromosome);

            final List<VariantContext> variantContexts = Lists.newArrayList();
            final Consumer<SageVariant> phasedConsumer = variant -> {
                if (shouldWrite(variant)) {
                    final VariantContext context = SageVariantContextFactory.create(variant);
                    variantContexts.add(context);
                }
            };

            final Phase phase = new Phase(reference, sageVariantFactory, phasedConsumer);
            for (CompletableFuture<SagePipelineData> region : regions) {
                for (SageVariant sageVariant : region.join().results()) {
                    phase.accept(sageVariant);
                }
            }

            phase.close(); //TODO: Rename to flush

            return variantContexts;
        });
    }

    private boolean shouldWrite(@NotNull final SageVariant entry) {
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
