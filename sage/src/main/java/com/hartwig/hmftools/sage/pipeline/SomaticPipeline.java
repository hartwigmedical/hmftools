package com.hartwig.hmftools.sage.pipeline;

import java.util.List;
import java.util.concurrent.CompletableFuture;
import java.util.concurrent.Executor;
import java.util.function.Supplier;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.genome.region.GenomeRegion;
import com.hartwig.hmftools.common.variant.hotspot.VariantHotspot;
import com.hartwig.hmftools.sage.config.SageConfig;
import com.hartwig.hmftools.sage.context.AltContext;
import com.hartwig.hmftools.sage.context.RefContext;
import com.hartwig.hmftools.sage.context.RefSequence;
import com.hartwig.hmftools.sage.evidence.NormalEvidence;
import com.hartwig.hmftools.sage.evidence.PrimaryEvidence;
import com.hartwig.hmftools.sage.sam.SamSlicerFactory;
import com.hartwig.hmftools.sage.variant.SageVariant;
import com.hartwig.hmftools.sage.variant.SageVariantFactory;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

import htsjdk.samtools.reference.IndexedFastaSequenceFile;

public class SomaticPipeline implements Supplier<CompletableFuture<List<SageVariant>>> {

    private static final Logger LOGGER = LogManager.getLogger(SomaticPipeline.class);

    private final GenomeRegion region;
    private final SageConfig config;
    private final Executor executor;
    private final RefSequence refSequence;
    private final SageVariantFactory variantFactory;
    private final SamSlicerFactory samSlicerFactory;
    private final List<VariantHotspot> hotspots;
    private final List<GenomeRegion> panelRegions;

    SomaticPipeline(final GenomeRegion region, final SageConfig config, final Executor executor, final IndexedFastaSequenceFile refGenome,
            final SamSlicerFactory samSlicerFactory, @NotNull final List<VariantHotspot> hotspots,
            @NotNull final List<GenomeRegion> panelRegions) {
        this.region = region;
        this.config = config;
        this.executor = executor;
        this.variantFactory = new SageVariantFactory(config.filter(), hotspots, panelRegions);
        this.samSlicerFactory = samSlicerFactory;
        this.refSequence = new RefSequence(region, refGenome);
        this.hotspots = hotspots;
        this.panelRegions = panelRegions;
    }

    @Override
    public CompletableFuture<List<SageVariant>> get() {
        return submit();
    }

    @NotNull
    public CompletableFuture<List<SageVariant>> submit() {

        final NormalEvidence normalEvidence = new NormalEvidence(config, samSlicerFactory);
        final PrimaryEvidence primaryEvidence = new PrimaryEvidence(config, hotspots, panelRegions, samSlicerFactory);

        final SomaticPipelineData somaticPipelineData = new SomaticPipelineData(config.reference(), config.tumor().size(), variantFactory);
        List<String> samples = config.tumor();
        List<String> bams = config.tumorBam();

        final List<CompletableFuture<List<AltContext>>> tumorFutures = Lists.newArrayList();
        for (int i = 0; i < samples.size(); i++) {
            final String sample = samples.get(i);
            final String bam = bams.get(i);

            CompletableFuture<List<AltContext>> candidateFuture =
                    CompletableFuture.supplyAsync(() -> primaryEvidence.get(sample, bam, refSequence, region), executor);

            tumorFutures.add(candidateFuture);
        }

        final CompletableFuture<Void> doneTumor = CompletableFuture.allOf(tumorFutures.toArray(new CompletableFuture[tumorFutures.size()]));

        final CompletableFuture<List<RefContext>> normalFuture = doneTumor.thenApply(aVoid -> {

            for (int i = 0; i < tumorFutures.size(); i++) {
                CompletableFuture<List<AltContext>> future = tumorFutures.get(i);
                somaticPipelineData.addTumor(i, future.join());
            }

            return normalEvidence.get(config.referenceBam(), refSequence, region, somaticPipelineData.normalCandidates());
        });

        return normalFuture.thenApply(aVoid -> {

            somaticPipelineData.addNormal(normalFuture.join());

            return somaticPipelineData.results();
        });
    }

}
