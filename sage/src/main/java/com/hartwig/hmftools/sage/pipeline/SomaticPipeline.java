package com.hartwig.hmftools.sage.pipeline;

import static java.util.concurrent.CompletableFuture.allOf;
import static java.util.concurrent.CompletableFuture.supplyAsync;

import java.util.List;
import java.util.concurrent.CompletableFuture;
import java.util.concurrent.Executor;
import java.util.function.Predicate;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.genome.region.GenomeRegion;
import com.hartwig.hmftools.common.variant.hotspot.VariantHotspot;
import com.hartwig.hmftools.sage.candidate.Candidate;
import com.hartwig.hmftools.sage.candidate.Candidates;
import com.hartwig.hmftools.sage.config.SageConfig;
import com.hartwig.hmftools.sage.evidence.CandidateEvidence;
import com.hartwig.hmftools.sage.evidence.ReadContextEvidence;
import com.hartwig.hmftools.sage.read.ReadContextCounter;
import com.hartwig.hmftools.sage.read.ReadContextCounters;
import com.hartwig.hmftools.sage.ref.RefSequence;
import com.hartwig.hmftools.sage.sam.SamSlicerFactory;
import com.hartwig.hmftools.sage.select.HotspotSelector;
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
            @NotNull final List<GenomeRegion> highConfidenceRegions) {
        this.config = config;
        this.executor = executor;
        final SamSlicerFactory samSlicerFactory = new SamSlicerFactory(config, panelRegions);
        this.hotspots = hotspots;
        this.panelRegions = panelRegions;
        this.highConfidenceRegions = highConfidenceRegions;
        this.candidateEvidence = new CandidateEvidence(config, hotspots, samSlicerFactory, refGenome);
        this.readContextEvidence = new ReadContextEvidence(config, samSlicerFactory, refGenome);
        this.refGenome = refGenome;
    }

    @NotNull
    public CompletableFuture<List<SageVariant>> variants(@NotNull final GenomeRegion region) {

        final CompletableFuture<RefSequence> refSequenceFuture = supplyAsync(() -> new RefSequence(region, refGenome), executor);

        // Scan tumors for set of initial candidates
        final CompletableFuture<List<Candidate>> doneCandidates = refSequenceFuture.thenCompose(refSequence -> {
            if (region.start() == 1) {
                LOGGER.info("Processing chromosome {}", region.chromosome());
            }
            LOGGER.debug("Processing initial candidates of {}:{}", region.chromosome(), region.start());

            final Candidates initialCandidates = new Candidates();
            final List<CompletableFuture<Void>> candidateFutures = Lists.newArrayList();

            for (int i = 0; i < config.tumor().size(); i++) {
                final String sample = config.tumor().get(i);
                final String sampleBam = config.tumorBam().get(i);

                final CompletableFuture<Void> candidateFuture =
                        refSequenceFuture.thenApply(x -> candidateEvidence.get(sample, sampleBam, refSequence, region))
                                .thenAccept(initialCandidates::add);

                candidateFutures.add(candidateFuture);
            }
            return allOf(candidateFutures.toArray(new CompletableFuture[candidateFutures.size()])).thenApply(y -> initialCandidates.candidates());
        });

        // Scan tumors for evidence
        final CompletableFuture<ReadContextCounters> doneTumor = doneCandidates.thenCompose(initialCandidates -> {
            LOGGER.debug("Scanning tumor for evidence in {}:{}", region.chromosome(), region.start());

            final ReadContextCounters result = new ReadContextCounters(config.primaryTumor(), initialCandidates);
            final List<CompletableFuture<Void>> tumorFutures = Lists.newArrayList();

            for (int i = 0; i < config.tumor().size(); i++) {
                final String sample = config.tumor().get(i);
                final String sampleBam = config.tumorBam().get(i);

                final CompletableFuture<Void> tumorFuture = CompletableFuture.completedFuture(this)
                        .thenApply(x -> readContextEvidence.get(region, initialCandidates, sample, sampleBam))
                        .thenAccept(result::addCounters);

                tumorFutures.add(tumorFuture);
            }

            return allOf(tumorFutures.toArray(new CompletableFuture[tumorFutures.size()])).thenApply(x -> result);
        });

        // Scan references for evidence
        final CompletableFuture<ReadContextCounters> doneNormal = doneTumor.thenCompose(tumorReadContextCounters -> {
            LOGGER.debug("Scanning reference for evidence in {}:{}", region.chromosome(), region.start());

            final Predicate<ReadContextCounter> hardFilter = hardFilterEvidence(hotspots);
            final List<Candidate> candidates = tumorReadContextCounters.candidates(hardFilter);
            final ReadContextCounters result = new ReadContextCounters(config.primaryReference(), candidates);

            final List<CompletableFuture<Void>> normalFutures = Lists.newArrayList();

            for (int i = 0; i < config.reference().size(); i++) {
                final String sample = config.reference().get(i);
                final String sampleBam = config.referenceBam().get(i);

                final CompletableFuture<Void> normalFuture = CompletableFuture.completedFuture(this)
                        .thenApply(x -> readContextEvidence.get(region, candidates, sample, sampleBam))
                        .thenAccept(result::addCounters);

                normalFutures.add(normalFuture);
            }

            return allOf(normalFutures.toArray(new CompletableFuture[normalFutures.size()])).thenApply(x -> result);
        });

        final CompletableFuture<List<SageVariant>> variantFuture =
                doneNormal.thenCombine(doneTumor, (normalCandidates, tumorCandidates) -> {
                    LOGGER.debug("Collating evidence in {}:{}", region.chromosome(), region.start());
                    final SageVariantFactory variantFactory =
                            new SageVariantFactory(config.filter(), hotspots, panelRegions, highConfidenceRegions);

                    final List<Candidate> candidates = normalCandidates.candidates(x -> true);

                    // Combine normal and tumor together and create variants
                    final List<SageVariant> result = Lists.newArrayList();
                    for (Candidate candidate : candidates) {
                        final List<ReadContextCounter> normal = normalCandidates.readContextCounters(candidate.variant());
                        final List<ReadContextCounter> tumor = tumorCandidates.readContextCounters(candidate.variant());
                        SageVariant sageVariant = variantFactory.createPairedTumorNormal(normal, tumor);
                        result.add(sageVariant);
                    }

                    return result;
                });

        return variantFuture;
    }

    @NotNull
    private Predicate<ReadContextCounter> hardFilterEvidence(@NotNull final List<VariantHotspot> variants) {
        return config.filter().hardFilter(new HotspotSelector(variants));
    }
}
