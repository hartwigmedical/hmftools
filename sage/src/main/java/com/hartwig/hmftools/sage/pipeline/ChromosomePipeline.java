package com.hartwig.hmftools.sage.pipeline;

import java.io.IOException;
import java.util.List;
import java.util.concurrent.CompletableFuture;
import java.util.function.Consumer;
import java.util.function.Function;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.genome.region.GenomeRegion;
import com.hartwig.hmftools.common.variant.hotspot.VariantHotspot;
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
    private final List<VariantHotspot> hotspots;
    private final List<GenomeRegion> panelRegions;

    ChromosomePipeline(@NotNull final String chromosome, @NotNull final SageConfig config,
            @NotNull final IndexedFastaSequenceFile reference, @NotNull final SageVariantFactory sageVariantFactory,
            @NotNull final Function<SageVariant, VariantContext> variantContextFactory, @NotNull final List<VariantHotspot> hotspots,
            @NotNull final List<GenomeRegion> panelRegions) throws IOException {
        this.chromosome = chromosome;
        this.config = config;
        this.reference = reference;
        this.sageVariantFactory = sageVariantFactory;
        this.variantContextFactory = variantContextFactory;
        this.sageVCF = new SageChromosomeVCF(chromosome, config);
        this.hotspots = hotspots;
        this.panelRegions = panelRegions;
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
        final Phase phase = new Phase(config, reference, sageVariantFactory, hotspots, panelRegions, phasedConsumer);

        final CompletableFuture<Void> done = CompletableFuture.allOf(regions.toArray(new CompletableFuture[regions.size()]));

        return done.thenApply(aVoid -> {
            LOGGER.info("Phasing chromosome {}", chromosome);

            // Phasing must be done in a separate thread as we re-query the MNVs
            for (final CompletableFuture<List<SageVariant>> region : regions) {
                region.join().forEach(phase);
            }

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
