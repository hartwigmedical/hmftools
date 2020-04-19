package com.hartwig.hmftools.sage.pipeline;

import static java.util.concurrent.CompletableFuture.supplyAsync;

import java.util.List;
import java.util.Map;
import java.util.concurrent.CompletableFuture;
import java.util.concurrent.Executor;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.genome.region.GenomeRegion;
import com.hartwig.hmftools.common.variant.hotspot.VariantHotspot;
import com.hartwig.hmftools.sage.candidate.Candidate;
import com.hartwig.hmftools.sage.candidate.Candidates;
import com.hartwig.hmftools.sage.config.SageConfig;
import com.hartwig.hmftools.sage.evidence.CandidateEvidence;
import com.hartwig.hmftools.sage.evidence.ReadContextEvidence;
import com.hartwig.hmftools.sage.quality.QualityRecalibrationMap;
import com.hartwig.hmftools.sage.read.ReadContextCounter;
import com.hartwig.hmftools.sage.read.ReadContextCounters;
import com.hartwig.hmftools.sage.ref.RefSequence;
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
    private final ReferenceSequenceFile refGenome;
    private final List<GenomeRegion> panelRegions;
    private final CandidateEvidence candidateEvidence;
    private final ReadContextEvidence readContextEvidence;
    private final List<GenomeRegion> highConfidenceRegions;

    SomaticPipeline(@NotNull final SageConfig config, @NotNull final Executor executor, @NotNull final ReferenceSequenceFile refGenome,
            @NotNull final List<VariantHotspot> hotspots, @NotNull final List<GenomeRegion> panelRegions,
            @NotNull final List<GenomeRegion> highConfidenceRegions,
            @NotNull final Map<String, QualityRecalibrationMap> qualityRecalibrationMap) {
        this.config = config;
        this.executor = executor;
        final SamSlicerFactory samSlicerFactory = new SamSlicerFactory(config, panelRegions);
        this.hotspots = hotspots;
        this.panelRegions = panelRegions;
        this.highConfidenceRegions = highConfidenceRegions;
        this.candidateEvidence = new CandidateEvidence(config, hotspots, panelRegions, samSlicerFactory, refGenome);
        this.readContextEvidence = new ReadContextEvidence(config, samSlicerFactory, refGenome, qualityRecalibrationMap);
        this.refGenome = refGenome;
    }

    @NotNull
    public CompletableFuture<List<SageVariant>> variants(@NotNull final GenomeRegion region) {

        final CompletableFuture<RefSequence> refSequenceFuture = supplyAsync(() -> new RefSequence(region, refGenome), executor);

        final CompletableFuture<List<Candidate>> initialCandidates = bamCandidates(region, refSequenceFuture);
        final CompletableFuture<ReadContextCounters> tumorEvidence = bamEvidence(region, config.tumor(), config.tumorBam(), initialCandidates);

        final CompletableFuture<List<Candidate>> finalCandidates = filteredCandidates(tumorEvidence);
        final CompletableFuture<ReadContextCounters> normalEvidence = bamEvidence(region, config.reference(), config.referenceBam(), finalCandidates);

        return combine(region, finalCandidates, tumorEvidence, normalEvidence);
    }

    @NotNull
    private CompletableFuture<List<Candidate>> bamCandidates(@NotNull final GenomeRegion region,
            final CompletableFuture<RefSequence> refSequenceFuture) {
        return refSequenceFuture.thenCompose(refSequence -> {
            if (region.start() == 1) {
                LOGGER.info("Processing chromosome {}", region.chromosome());
            }
            LOGGER.debug("Processing candidates in {}:{}", region.chromosome(), region.start());

            final Candidates initialCandidates = new Candidates(hotspots, panelRegions, highConfidenceRegions);

            CompletableFuture<Void> done = CompletableFuture.completedFuture(null);
            for (int i = 0; i < config.tumor().size(); i++) {
                final String sample = config.tumor().get(i);
                final String sampleBam = config.tumorBam().get(i);
                done = done.thenApply(aVoid -> candidateEvidence.get(sample, sampleBam, refSequence, region))
                        .thenAccept(initialCandidates::add);
            }
            return done.thenApply(y -> initialCandidates.candidates());
        });
    }

    @NotNull
    private CompletableFuture<ReadContextCounters> bamEvidence(@NotNull final GenomeRegion region, @NotNull final List<String> samples,
            @NotNull final List<String> sampleBams, @NotNull final CompletableFuture<List<Candidate>> candidates) {
        // Scan tumors for evidence
        return candidates.thenCompose(initialCandidates -> {
            LOGGER.debug("Scanning for evidence in {}:{}", region.chromosome(), region.start());
            final String primarySample = samples.isEmpty() ? "PRIMARY" : samples.get(0);

            final ReadContextCounters result = new ReadContextCounters(primarySample, initialCandidates);

            CompletableFuture<Void> done = CompletableFuture.completedFuture(null);
            for (int i = 0; i < samples.size(); i++) {
                final String sample = samples.get(i);
                final String sampleBam = sampleBams.get(i);

                done = done.thenApply(x -> readContextEvidence.get(initialCandidates, sample, sampleBam)).thenAccept(result::addCounters);
            }

            return done.thenApply(x -> result);
        });
    }

    @NotNull
    private CompletableFuture<List<SageVariant>> combine(@NotNull final GenomeRegion region,
            final CompletableFuture<List<Candidate>> candidates, final CompletableFuture<ReadContextCounters> doneTumor,
            final CompletableFuture<ReadContextCounters> doneNormal) {
        return doneNormal.thenCombine(doneTumor, (normalCandidates, tumorCandidates) -> {
            LOGGER.info("Gathering evidence in {}:{}", region.chromosome(), region.start());
            final SageVariantFactory variantFactory = new SageVariantFactory(config.filter());

            // Combine normal and tumor together and create variants
            final List<SageVariant> result = Lists.newArrayList();
            for (Candidate candidate : candidates.join()) {
                final List<ReadContextCounter> normal = normalCandidates.readContextCounters(candidate.variant());
                final List<ReadContextCounter> tumor = tumorCandidates.readContextCounters(candidate.variant());
                SageVariant sageVariant = variantFactory.create(normal, tumor);
                result.add(sageVariant);
            }

            return result;
        });
    }

    @NotNull
    private CompletableFuture<List<Candidate>> filteredCandidates(final CompletableFuture<ReadContextCounters> tumorEvidence) {
        return tumorEvidence.thenApply(x -> x.candidates(config.filter().readContextFilter()));
    }
}
