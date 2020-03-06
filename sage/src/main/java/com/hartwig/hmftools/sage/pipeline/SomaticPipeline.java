package com.hartwig.hmftools.sage.pipeline;

import java.util.List;
import java.util.concurrent.CompletableFuture;
import java.util.concurrent.Executor;

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

import htsjdk.samtools.reference.ReferenceSequenceFile;

public class SomaticPipeline implements SageVariantPipeline {

    private static final Logger LOGGER = LogManager.getLogger(SomaticPipeline.class);

    private final SageConfig config;
    private final Executor executor;
    private final List<VariantHotspot> hotspots;
    private final List<GenomeRegion> panelRegions;
    private final List<GenomeRegion> highConfidenceRegions;
    private final PrimaryEvidence primaryEvidence;
    private final NormalEvidence normalEvidence;
    private final ReferenceSequenceFile refGenome;

    SomaticPipeline(@NotNull final SageConfig config, @NotNull final Executor executor, @NotNull final ReferenceSequenceFile refGenome,
            @NotNull final List<VariantHotspot> hotspots, @NotNull final List<GenomeRegion> panelRegions,
            @NotNull final List<GenomeRegion> highConfidenceRegions) {
        this.config = config;
        this.executor = executor;
        final SamSlicerFactory samSlicerFactory = new SamSlicerFactory(config, panelRegions);
        this.hotspots = hotspots;
        this.panelRegions = panelRegions;
        this.highConfidenceRegions = highConfidenceRegions;
        this.primaryEvidence = new PrimaryEvidence(config, hotspots, samSlicerFactory, refGenome);
        this.normalEvidence = new NormalEvidence(config, samSlicerFactory, refGenome);
        this.refGenome = refGenome;
    }

    @NotNull
    public CompletableFuture<List<SageVariant>> variants(@NotNull final GenomeRegion region) {

        final SageVariantFactory variantFactory = new SageVariantFactory(config.filter(), hotspots, panelRegions, highConfidenceRegions);
        final SomaticPipelineData somaticPipelineData = new SomaticPipelineData(config, variantFactory);
        List<String> samples = config.tumor();
        List<String> bams = config.tumorBam();

        final CompletableFuture<RefSequence> refSequenceFuture =
                CompletableFuture.supplyAsync(() -> new RefSequence(region, refGenome), executor);

        final List<CompletableFuture<List<AltContext>>> tumorFutures = Lists.newArrayList();
        final CompletableFuture<Void> doneTumor = refSequenceFuture.thenCompose(refSequence -> {
            for (int i = 0; i < samples.size(); i++) {
                final String sample = samples.get(i);
                final String bam = bams.get(i);

                CompletableFuture<List<AltContext>> candidateFuture =
                        CompletableFuture.supplyAsync(() -> primaryEvidence.get(sample, bam, refSequence, region), executor);

                tumorFutures.add(candidateFuture);
            }

            return CompletableFuture.allOf(tumorFutures.toArray(new CompletableFuture[tumorFutures.size()]));
        });

        final CompletableFuture<Void> doneMergedTumor = doneTumor.thenCompose(aVoid -> {
            for (int i = 0; i < tumorFutures.size(); i++) {
                CompletableFuture<List<AltContext>> future = tumorFutures.get(i);
                somaticPipelineData.addTumor(i, future.join());
            }
            return CompletableFuture.completedFuture(null);
        });

        final List<CompletableFuture<List<RefContext>>> normalFutures = Lists.newArrayList();
        final CompletableFuture<Void> doneNormal = doneMergedTumor.thenCompose(aVoid -> {
            for (int i = 0; i < config.reference().size(); i++) {
                final String sample = config.reference().get(i);
                final String sampleBam = config.referenceBam().get(i);

                CompletableFuture<List<RefContext>> normalFuture =
                        CompletableFuture.supplyAsync(() -> normalEvidence.get(refSequenceFuture.join(),
                                region,
                                somaticPipelineData.normalCandidates(sample),
                                sampleBam));

                normalFutures.add(normalFuture);
            }
            return CompletableFuture.allOf(normalFutures.toArray(new CompletableFuture[normalFutures.size()]));
        });

        final CompletableFuture<Void> doneMergedNormal = doneNormal.thenCompose(aVoid -> {
            for (int i = 0; i < normalFutures.size(); i++) {
                CompletableFuture<List<RefContext>> future = normalFutures.get(i);
                somaticPipelineData.addNormal(i, future.join());
            }
            return CompletableFuture.completedFuture(null);
        });

        return doneMergedNormal.thenApply(aVoid -> somaticPipelineData.results());
    }
}
