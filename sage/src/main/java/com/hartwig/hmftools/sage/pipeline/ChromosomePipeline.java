package com.hartwig.hmftools.sage.pipeline;

import java.io.IOException;
import java.util.List;
import java.util.concurrent.CompletableFuture;
import java.util.function.Consumer;
import java.util.function.Function;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.sage.config.SageConfig;
import com.hartwig.hmftools.sage.context.AltContext;
import com.hartwig.hmftools.sage.phase.Phase;
import com.hartwig.hmftools.sage.variant.SageVariant;
import com.hartwig.hmftools.sage.variant.SageVariantFactory;
import com.hartwig.hmftools.sage.vcf.SageChromosomeVCF;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import htsjdk.variant.variantcontext.VariantContext;

public class ChromosomePipeline implements Consumer<CompletableFuture<List<SageVariant>>> {

    private static final Logger LOGGER = LogManager.getLogger(ChromosomePipeline.class);

    private final String chromosome;
    private final SageConfig config;
    private final SageChromosomeVCF sageVCF;
    private final IndexedFastaSequenceFile reference;
    private final SageVariantFactory sageVariantFactory;
    private final Function<SageVariant, VariantContext> variantContextFactory;
    private final List<CompletableFuture<List<SageVariant>>> regions = Lists.newArrayList();

    ChromosomePipeline(@NotNull final String chromosome, @NotNull final SageConfig config,
            @NotNull final IndexedFastaSequenceFile reference, @NotNull final SageVariantFactory sageVariantFactory,
            @NotNull final Function<SageVariant, VariantContext> variantContextFactory) throws IOException {
        this.chromosome = chromosome;
        this.config = config;
        this.reference = reference;
        this.sageVariantFactory = sageVariantFactory;
        this.variantContextFactory = variantContextFactory;
        this.sageVCF = new SageChromosomeVCF(chromosome, config);
    }

    @NotNull
    public String chromosome() {
        return chromosome;
    }

    @NotNull
    public String vcfFilename() {
        return sageVCF.filename();
    }

    @Override
    public void accept(final CompletableFuture<List<SageVariant>> regionData) {
        regions.add(regionData);
    }

    @NotNull
    public CompletableFuture<ChromosomePipeline> submit() {

        final Consumer<SageVariant> phasedConsumer = variant -> {
            if (include(variant)) {
                final VariantContext context = variantContextFactory.apply(variant);
                sageVCF.write(context);
            }
        };
        final Phase phase = new Phase(config, reference, sageVariantFactory, phasedConsumer);

        // Phasing must be done in (positional) order but we can do it eagerly as each new region comes in.
        // It is not necessary to wait for the entire chromosome to be finished to start.
        CompletableFuture<Void> done = CompletableFuture.completedFuture(null);
        for (final CompletableFuture<List<SageVariant>> region : regions) {
            done = done.thenCombine(region, (aVoid, sageVariants) -> {
                sageVariants.forEach(phase);
                return null;
            });
        }

        return done.thenApply(aVoid -> {
            phase.flush();
            sageVCF.close();
            LOGGER.info("Finished processing chromosome {}", chromosome);
            return ChromosomePipeline.this;
        });

    }

    private boolean include(@NotNull final SageVariant entry) {
        if (entry.isPassing()) {
            return true;
        }

        if (config.filter().hardFilter()) {
            return false;
        }

        if (config.germlineOnly()) {
            return true;
        }

        final AltContext normal = entry.normal();
        if (normal.rawAltSupport() > config.filter().hardMaxNormalAltSupport()) {
            return false;
        }

        return entry.primaryTumor().primaryReadContext().tumorQuality() >= config.filter().hardMinTumorQualFiltered();
    }

}
