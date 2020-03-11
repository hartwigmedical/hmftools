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
import com.hartwig.hmftools.sage.candidate.AltContextCandidates;
import com.hartwig.hmftools.sage.candidate.ReadContextCandidates;
import com.hartwig.hmftools.sage.config.SageConfig;
import com.hartwig.hmftools.sage.context.AltContext;
import com.hartwig.hmftools.sage.context.FixedRefContextCandidates;
import com.hartwig.hmftools.sage.context.FixedRefContextCandidatesFactory;
import com.hartwig.hmftools.sage.context.RefSequence;
import com.hartwig.hmftools.sage.evidence.CandidateEvidence;
import com.hartwig.hmftools.sage.evidence.StandardEvidence;
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
    private final List<GenomeRegion> panelRegions;
    private final List<GenomeRegion> highConfidenceRegions;
    private final CandidateEvidence initialEvidence;
    private final StandardEvidence normalEvidence;
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
        this.initialEvidence = new CandidateEvidence(config, hotspots, samSlicerFactory, refGenome);
        this.normalEvidence = new StandardEvidence(config, samSlicerFactory, refGenome);
        this.refGenome = refGenome;
    }

    @NotNull
    public CompletableFuture<List<SageVariant>> variants(@NotNull final GenomeRegion region) {

        final SageVariantFactory variantFactory = new SageVariantFactory(config.filter(), hotspots, panelRegions, highConfidenceRegions);

        final CompletableFuture<RefSequence> refSequenceFuture = supplyAsync(() -> new RefSequence(region, refGenome), executor);

        // Scan tumors for set of initial candidates
        final CompletableFuture<ReadContextCandidates> doneCandidates = refSequenceFuture.thenCompose(refSequence -> {

            if (region.start() == 1) {
                LOGGER.info("Processing chromosome {}", region.chromosome());
            }

            final List<CompletableFuture<Void>> candidateFutures = Lists.newArrayList();
            final ReadContextCandidates readContextCandidates = new ReadContextCandidates();

            for (int i = 0; i < config.tumor().size(); i++) {
                final String sample = config.tumor().get(i);
                final String sampleBam = config.tumorBam().get(i);

                final CompletableFuture<Void> candidateFuture =
                        supplyAsync(() -> initialEvidence.get(sample, sampleBam, refSequence, region), executor).thenAccept(
                                readContextCandidates::addCandidates);

                candidateFutures.add(candidateFuture);
            }
            return allOf(candidateFutures.toArray(new CompletableFuture[candidateFutures.size()])).thenApply(x -> readContextCandidates);
        });

        // Scan tumors for evidence
        final CompletableFuture<AltContextCandidates> doneTumor = doneCandidates.thenCompose(readContextCandidates -> {
            final FixedRefContextCandidatesFactory candidatesFactory = readContextCandidates.candidateFactory();
            final List<CompletableFuture<Void>> tumorFutures = Lists.newArrayList();
            final AltContextCandidates tumorCandidates = new AltContextCandidates(config.primaryTumor(), candidatesFactory.loci());

            for (int i = 0; i < config.tumor().size(); i++) {
                final String sample = config.tumor().get(i);
                final String sampleBam = config.tumorBam().get(i);
                final FixedRefContextCandidates fixedCandidates = candidatesFactory.create(sample);

                final CompletableFuture<Void> tumorFuture =
                        supplyAsync(() -> normalEvidence.get(refSequenceFuture.join(), region, fixedCandidates, sampleBam)).thenAccept(
                                tumorCandidates::addRefCandidate);

                tumorFutures.add(tumorFuture);
            }

            return allOf(tumorFutures.toArray(new CompletableFuture[tumorFutures.size()])).thenApply(x -> tumorCandidates);
        });

        // Scan references for evidence
        final CompletableFuture<AltContextCandidates> doneNormal = doneTumor.thenCompose(tumorCandidates -> {
            final Predicate<AltContext> hardFilter = hardFilterEvidence(hotspots);
            final FixedRefContextCandidatesFactory candidatesFactory = tumorCandidates.createFactory(hardFilter);
            final List<CompletableFuture<Void>> normalFutures = Lists.newArrayList();
            final AltContextCandidates normalCandidates = new AltContextCandidates(config.primaryReference(), candidatesFactory.loci());

            for (int i = 0; i < config.reference().size(); i++) {
                final String sample = config.reference().get(i);
                final String sampleBam = config.referenceBam().get(i);
                final FixedRefContextCandidates fixedCandidates = candidatesFactory.create(sample);

                final CompletableFuture<Void> normalFuture =
                        supplyAsync(() -> normalEvidence.get(refSequenceFuture.join(), region, fixedCandidates, sampleBam)).thenAccept(
                                normalCandidates::addRefCandidate);

                normalFutures.add(normalFuture);
            }

            return allOf(normalFutures.toArray(new CompletableFuture[normalFutures.size()])).
                    thenApply(x -> normalCandidates);
        });

        return doneNormal.thenApply(normalCandidates -> {

            // Combine normal and tumor together and create variants
            final List<SageVariant> result = Lists.newArrayList();
            final FixedRefContextCandidatesFactory sorted = normalCandidates.createFactory(x -> true);
            for (VariantHotspot variant : sorted.loci()) {
                final List<AltContext> normal = normalCandidates.altContexts(variant);
                final List<AltContext> tumor = doneTumor.join().altContexts(variant);

                SageVariant sageVariant = variantFactory.create(normal, tumor);
                result.add(sageVariant);
            }

            return result;
        });
    }

    @NotNull
    private Predicate<AltContext> hardFilterEvidence(@NotNull final List<VariantHotspot> hotspots) {
        final HotspotSelector tierSelector = new HotspotSelector(hotspots);
        return altContext -> altContext.primaryReadContext().tumorQuality() >= config.filter().hardMinTumorQual() || tierSelector.isHotspot(
                altContext);
    }
}
