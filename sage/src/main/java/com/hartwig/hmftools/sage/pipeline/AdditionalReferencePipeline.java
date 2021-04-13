package com.hartwig.hmftools.sage.pipeline;

import static java.util.concurrent.CompletableFuture.supplyAsync;

import java.util.List;
import java.util.Map;
import java.util.concurrent.CompletableFuture;
import java.util.concurrent.Executor;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.genome.region.GenomeRegion;
import com.hartwig.hmftools.common.variant.hotspot.VariantHotspot;
import com.hartwig.hmftools.sage.candidate.Candidate;
import com.hartwig.hmftools.sage.candidate.CandidateSerialization;
import com.hartwig.hmftools.sage.config.SageConfig;
import com.hartwig.hmftools.sage.quality.QualityRecalibrationMap;
import com.hartwig.hmftools.sage.read.ReadContextCounter;
import com.hartwig.hmftools.sage.read.ReadContextCounters;
import com.hartwig.hmftools.sage.ref.RefSequence;
import com.hartwig.hmftools.sage.variant.SageVariantContextFactory;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

import htsjdk.samtools.reference.ReferenceSequenceFile;
import htsjdk.variant.variantcontext.VariantContext;

public class AdditionalReferencePipeline {

    private static final Logger LOGGER = LogManager.getLogger(AdditionalReferencePipeline.class);

    private final SageConfig config;
    private final EvidenceStage evidenceStage;
    private final ReferenceSequenceFile refGenome;
    private final Executor executor;

    public AdditionalReferencePipeline(@NotNull final SageConfig config, @NotNull final Executor executor, ReferenceSequenceFile refGenome,
            @NotNull final Map<String, QualityRecalibrationMap> qualityRecalibrationMap) {
        this.config = config;
        this.refGenome = refGenome;
        this.evidenceStage = new EvidenceStage(config, refGenome, qualityRecalibrationMap);
        this.executor = executor;
    }

    @NotNull
    public CompletableFuture<List<VariantContext>> appendReference(@NotNull final GenomeRegion region,
            @NotNull final List<VariantContext> variants) {
        if (variants.isEmpty()) {
            return CompletableFuture.completedFuture(Lists.newArrayList());
        }

        final CompletableFuture<RefSequence> refSequenceFuture = supplyAsync(() -> {
            if (region.start() == 1) {
                LOGGER.info("Processing chromosome {}", region.chromosome());
            }
            return new RefSequence(region, refGenome);
        }, executor);

        final CompletableFuture<List<Candidate>> candidateFutures = refSequenceFuture.thenApply(x -> variants.stream()
                .map(y -> CandidateSerialization.toCandidate(y, x))
                .collect(Collectors.toList()));

        final CompletableFuture<ReadContextCounters> evidenceFutures =
                evidenceStage.evidence(config.reference(), config.referenceBam(), candidateFutures);

        return evidenceFutures.thenApply(x -> update(x, variants));
    }

    public List<VariantContext> update(final ReadContextCounters readContextCounters, final List<VariantContext> variantContexts) {
        final List<VariantContext> result = Lists.newArrayList();
        for (VariantContext old : variantContexts) {
            VariantHotspot variant = CandidateSerialization.toVariantHotspot(old);
            List<ReadContextCounter> counters = readContextCounters.readContextCounters(variant);
            result.add(SageVariantContextFactory.addGenotype(old, counters));
        }

        return result;
    }
}
