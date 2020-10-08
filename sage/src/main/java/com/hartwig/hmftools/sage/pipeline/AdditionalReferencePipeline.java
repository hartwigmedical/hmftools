package com.hartwig.hmftools.sage.pipeline;

import static java.util.concurrent.CompletableFuture.supplyAsync;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
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
import com.hartwig.hmftools.sage.variant.SageVariantContextFactory;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import htsjdk.variant.variantcontext.VariantContext;

public class AdditionalReferencePipeline implements AutoCloseable {

    private static final Logger LOGGER = LogManager.getLogger(AdditionalReferencePipeline.class);

    private final String contig;
    private final SageConfig config;
    private final EvidenceStage evidenceStage;
    private final IndexedFastaSequenceFile refGenome;
    private final Executor executor;
    private final CandidateSerialization serialization;

    public AdditionalReferencePipeline(@NotNull final String contig, @NotNull final SageConfig config, @NotNull final Executor executor,
            @NotNull final Map<String, QualityRecalibrationMap> qualityRecalibrationMap) throws FileNotFoundException {
        this.contig = contig;
        this.config = config;
        this.refGenome = new IndexedFastaSequenceFile(new File(config.refGenome()));
        this.evidenceStage = new EvidenceStage(config, refGenome, qualityRecalibrationMap);
        this.executor = executor;
        this.serialization = new CandidateSerialization(refGenome);
    }

    @NotNull
    public CompletableFuture<List<VariantContext>> appendReference(@NotNull final GenomeRegion region,
            @NotNull final List<VariantContext> variants) {


        final CompletableFuture<List<Candidate>> candidateFutures =
                CompletableFuture.completedFuture(variants.stream().map(serialization::toCandidate).collect(Collectors.toList()));

        final CompletableFuture<List<VariantContext>> existingFuture = supplyAsync(() -> {
            if (region.start() == 1) {
                LOGGER.info("Processing chromosome {}", contig);
            }
            return variants;
        }, executor);

        final CompletableFuture<ReadContextCounters> evidenceFutures =
                evidenceStage.evidence(config.reference(), config.referenceBam(), candidateFutures);

        return evidenceFutures.thenCombine(existingFuture, this::update);

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

    @Override
    public void close() throws IOException {
        refGenome.close();
    }
}
