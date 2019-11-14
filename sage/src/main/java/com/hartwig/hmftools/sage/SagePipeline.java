package com.hartwig.hmftools.sage;

import java.util.List;
import java.util.concurrent.CompletableFuture;
import java.util.concurrent.Executor;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.genome.region.GenomeRegion;
import com.hartwig.hmftools.sage.config.SageConfig;
import com.hartwig.hmftools.sage.context.AltContext;
import com.hartwig.hmftools.sage.context.NormalRefContextSupplier;
import com.hartwig.hmftools.sage.context.RefContext;
import com.hartwig.hmftools.sage.context.RefSequence;
import com.hartwig.hmftools.sage.context.TumorAltContextSupplier;
import com.hartwig.hmftools.sage.sam.SamSlicerFactory;
import com.hartwig.hmftools.sage.variant.SageVariant;
import com.hartwig.hmftools.sage.variant.SageVariantFactory;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

import htsjdk.samtools.reference.IndexedFastaSequenceFile;

public class SagePipeline {

    private static final Logger LOGGER = LogManager.getLogger(SagePipeline.class);

    private final GenomeRegion region;
    private final SageConfig config;
    private final Executor executor;
    private final RefSequence refSequence;
    private final SageVariantFactory variantFactory;
    private final SamSlicerFactory samSlicerFactory;

    public SagePipeline(final GenomeRegion region, final SageConfig config, final Executor executor,
            final IndexedFastaSequenceFile refGenome, final SageVariantFactory variantFactory, final SamSlicerFactory samSlicerFactory) {
        this.region = region;
        this.config = config;
        this.executor = executor;
        this.variantFactory = variantFactory;
        this.samSlicerFactory = samSlicerFactory;
        this.refSequence = new RefSequence(region, refGenome);
    }

    @NotNull
    public CompletableFuture<List<SageVariant>> submit() {

        final SagePipelineData sagePipelineData = new SagePipelineData(config.reference(), config.tumor().size(), variantFactory);
        List<String> samples = config.tumor();
        List<String> bams = config.tumorBam();

        final List<CompletableFuture<List<AltContext>>> tumorFutures = Lists.newArrayList();
        for (int i = 0; i < samples.size(); i++) {
            final String sample = samples.get(i);
            final String bam = bams.get(i);

            CompletableFuture<List<AltContext>> candidateFuture =
                    CompletableFuture.supplyAsync(new TumorAltContextSupplier(config, sample, region, bam, refSequence, samSlicerFactory), executor);

            tumorFutures.add(candidateFuture);
        }

        final CompletableFuture<Void> doneTumor = CompletableFuture.allOf(tumorFutures.toArray(new CompletableFuture[tumorFutures.size()]));

        final CompletableFuture<List<RefContext>> normalFuture = doneTumor.thenApply(aVoid -> {

            for (int i = 0; i < tumorFutures.size(); i++) {
                CompletableFuture<List<AltContext>> future = tumorFutures.get(i);
                sagePipelineData.addTumor(i, future.join());
            }

            return new NormalRefContextSupplier(config,
                    region,
                    config.referenceBam(),
                    refSequence,
                    sagePipelineData.normalCandidates(), samSlicerFactory).get();
        });

        return normalFuture.thenApply(aVoid -> {

            sagePipelineData.addNormal(normalFuture.join());

            return sagePipelineData.results();
        });
    }
}
