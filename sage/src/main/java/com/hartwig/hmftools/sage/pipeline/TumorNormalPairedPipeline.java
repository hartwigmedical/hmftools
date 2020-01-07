package com.hartwig.hmftools.sage.pipeline;

import java.io.IOException;
import java.util.List;
import java.util.concurrent.CompletableFuture;
import java.util.concurrent.Executor;
import java.util.function.Function;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.genome.region.GenomeRegion;
import com.hartwig.hmftools.common.variant.hotspot.VariantHotspot;
import com.hartwig.hmftools.sage.config.SageConfig;
import com.hartwig.hmftools.sage.context.AltContext;
import com.hartwig.hmftools.sage.variant.SageVariant;
import com.hartwig.hmftools.sage.variant.SageVariantContextFactory;
import com.hartwig.hmftools.sage.variant.SageVariantFactory;
import com.hartwig.hmftools.sage.vcf.SageChromosomeVCF;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import htsjdk.variant.variantcontext.VariantContext;

public class TumorNormalPairedPipeline {

    private static final Logger LOGGER = LogManager.getLogger(TumorNormalPipeline.class);


    private final String chromosome;
    private final SageConfig config;
    private final IndexedFastaSequenceFile refGenome;
    private final SageChromosomeVCF sageVCF;
    private final Executor executor;
    private final List<VariantHotspot> hotspots;
    private final List<GenomeRegion> panelRegions;
    private final SageVariantFactory variantFactory;
    private final Function<SageVariant, VariantContext> variantContextFactory;

    public TumorNormalPairedPipeline(final String chromosome, final SageConfig config, final IndexedFastaSequenceFile refGenome,
            @NotNull final List<VariantHotspot> hotspots, @NotNull final List<GenomeRegion> panelRegions, final Executor executor)
            throws IOException {
        this.chromosome = chromosome;
        this.config = config;
        this.sageVCF = new SageChromosomeVCF(chromosome, config);
        this.executor = executor;
        this.hotspots = hotspots;
        this.panelRegions = panelRegions;
        this.refGenome = refGenome;
        this.variantFactory = new SageVariantFactory(config.filter(), hotspots, panelRegions);
        this.variantContextFactory = SageVariantContextFactory::pairedTumorNormal;

    }

    @NotNull
    public String chromosome() {
        return chromosome;
    }

    @NotNull
    public String vcfFilename() {
        return sageVCF.filename();
    }

    @NotNull
    public CompletableFuture<TumorNormalPairedPipeline> submit() {
        return submit(1, refGenome.getSequence(chromosome).length());
    }

    @NotNull
    public CompletableFuture<TumorNormalPairedPipeline> submit(int minPosition, int maxPosition) {

        final SomaticPipelineData somaticPipelineData = new SomaticPipelineData(config.reference(), config.tumor().size(), variantFactory);

        final List<CompletableFuture<TumorPipeline>> tumorFutures = Lists.newArrayList();
        for (int i = 0; i < config.tumor().size(); i++) {
            tumorFutures.add(new TumorPipeline(0, chromosome, config, refGenome, hotspots, panelRegions, executor).submit(minPosition,
                    maxPosition));
        }

        final CompletableFuture<Void> doneTumor = CompletableFuture.allOf(tumorFutures.toArray(new CompletableFuture[tumorFutures.size()]));

        final CompletableFuture<TumorNormalPipeline> normalFuture = doneTumor.thenCompose(aVoid -> {

            for (int i = 0; i < tumorFutures.size(); i++) {
                CompletableFuture<TumorPipeline> future = tumorFutures.get(i);
                somaticPipelineData.addTumor(i, future.join().candidates());
            }

            final TumorNormalPipeline tumorNormalPipeline =
                    new TumorNormalPipeline(chromosome, config, refGenome, panelRegions, executor, somaticPipelineData.normalCandidates());
            return tumorNormalPipeline.submit(minPosition, maxPosition);
        });

        return normalFuture.thenApply(aVoid -> {

            LOGGER.info("Combining1");

            somaticPipelineData.addNormal(normalFuture.join().candidates());
            LOGGER.info("Combining2");
            somaticPipelineData.results(variantFactory)
                    .stream()
                    .filter(this::include)
                    .peek(x -> x.localPhaseSet(x.primaryTumor().localPhaseSet()))
                    .map(variantContextFactory)
                    .forEach(sageVCF::write);

            sageVCF.close();
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
