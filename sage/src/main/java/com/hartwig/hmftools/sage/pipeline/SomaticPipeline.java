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
import com.hartwig.hmftools.sage.candidate.AltContextCollection;
import com.hartwig.hmftools.sage.candidate.ReadContextCollection;
import com.hartwig.hmftools.sage.config.SageConfig;
import com.hartwig.hmftools.sage.context.AltContext;
import com.hartwig.hmftools.sage.context.RefContextFixedFactory;
import com.hartwig.hmftools.sage.context.RefContextFixedFactoryForSample;
import com.hartwig.hmftools.sage.context.RefSequence;
import com.hartwig.hmftools.sage.evidence.CandidateEvidence;
import com.hartwig.hmftools.sage.evidence.FixedEvidence;
import com.hartwig.hmftools.sage.read.ReadContextCounter;
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
    private final FixedEvidence normalEvidence;
    private final CandidateEvidence candidateEvidence;
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
        this.normalEvidence = new FixedEvidence(config, samSlicerFactory, refGenome);
        this.refGenome = refGenome;
    }

    @NotNull
    public CompletableFuture<List<SageVariant>> variants(@NotNull final GenomeRegion region) {

        final CompletableFuture<RefSequence> refSequenceFuture = supplyAsync(() -> new RefSequence(region, refGenome), executor);

        // Scan tumors for set of initial candidates
        final CompletableFuture<ReadContextCollection> doneCandidates = refSequenceFuture.thenCompose(refSequence -> {
            if (region.start() == 1) {
                LOGGER.info("Processing chromosome {}", region.chromosome());
            }
            LOGGER.debug("Processing initial candidates of {}:{}", region.chromosome(), region.start());

            final List<CompletableFuture<Void>> candidateFutures = Lists.newArrayList();
            final ReadContextCollection readContextCandidates = new ReadContextCollection(config);

            for (int i = 0; i < config.tumor().size(); i++) {
                final String sample = config.tumor().get(i);
                final String sampleBam = config.tumorBam().get(i);

                final CompletableFuture<Void> candidateFuture =
                        refSequenceFuture.thenApply(x -> candidateEvidence.get(sample, sampleBam, refSequence, region))
                                .thenAccept(readContextCandidates::addCandidates);

                candidateFutures.add(candidateFuture);
            }
            return allOf(candidateFutures.toArray(new CompletableFuture[candidateFutures.size()])).thenApply(y -> readContextCandidates);
        });

        // Scan tumors for evidence
        final CompletableFuture<AltContextCollection> doneTumor = doneCandidates.thenCompose(readContextCandidates -> {
            LOGGER.info("Scanning tumor for evidence in {}:{}", region.chromosome(), region.start());

            final RefContextFixedFactoryForSample candidatesFactory = readContextCandidates.supplier();
            final List<CompletableFuture<Void>> tumorFutures = Lists.newArrayList();
            final AltContextCollection tumorCandidates = new AltContextCollection(config, config.primaryTumor(), candidatesFactory.loci());

            for (int i = 0; i < config.tumor().size(); i++) {
                final String sample = config.tumor().get(i);
                final String sampleBam = config.tumorBam().get(i);
                final RefContextFixedFactory fixedCandidates = candidatesFactory.create(sample);

                final CompletableFuture<Void> tumorFuture = CompletableFuture.completedFuture(this)
                        .thenApply(x -> normalEvidence.get(region, fixedCandidates, sampleBam))
                        .thenAccept(tumorCandidates::addRefContexts);

                tumorFutures.add(tumorFuture);
            }

            return allOf(tumorFutures.toArray(new CompletableFuture[tumorFutures.size()])).thenApply(x -> tumorCandidates);
        });

        // Scan references for evidence
        final CompletableFuture<AltContextCollection> doneNormal = doneTumor.thenCompose(tumorCandidates -> {
            LOGGER.debug("Scanning reference for evidence in {}:{}", region.chromosome(), region.start());

            final Predicate<AltContext> hardFilter = hardFilterEvidence(hotspots);
            final RefContextFixedFactoryForSample candidatesFactory = tumorCandidates.createFactory(hardFilter);
            final List<CompletableFuture<Void>> normalFutures = Lists.newArrayList();
            final AltContextCollection normalCandidates =
                    new AltContextCollection(config, config.primaryReference(), candidatesFactory.loci());

            for (int i = 0; i < config.reference().size(); i++) {
                final String sample = config.reference().get(i);
                final String sampleBam = config.referenceBam().get(i);
                final RefContextFixedFactory fixedCandidates = candidatesFactory.create(sample);

                final CompletableFuture<Void> normalFuture = CompletableFuture.completedFuture(this)
                        .thenApply(x -> normalEvidence.get(region, fixedCandidates, sampleBam))
                        .thenAccept(normalCandidates::addRefContexts);

                normalFutures.add(normalFuture);
            }

            return allOf(normalFutures.toArray(new CompletableFuture[normalFutures.size()])).
                    thenApply(x -> normalCandidates);
        });

        final CompletableFuture<List<SageVariant>> variantFuture =
                doneNormal.thenCombine(doneTumor, (normalCandidates, tumorCandidates) -> {
                    LOGGER.debug("Collating evidence in {}:{}", region.chromosome(), region.start());
                    final SageVariantFactory variantFactory =
                            new SageVariantFactory(config.filter(), hotspots, panelRegions, highConfidenceRegions);

                    // Combine normal and tumor together and create variants
                    final List<SageVariant> result = Lists.newArrayList();
                    final RefContextFixedFactoryForSample sorted = normalCandidates.createFactory(x -> true);
                    for (VariantHotspot variant : sorted.loci()) {
                        final List<ReadContextCounter> normal = normalCandidates.altContexts(variant);
                        final List<ReadContextCounter> tumor = tumorCandidates.altContexts(variant);

                        SageVariant sageVariant = variantFactory.create(normal, tumor);
                        result.add(sageVariant);
                    }

                    return result;
                });

        return variantFuture;
    }

    @NotNull
    private Predicate<AltContext> hardFilterEvidence(@NotNull final List<VariantHotspot> hotspots) {
        final HotspotSelector tierSelector = new HotspotSelector(hotspots);
        return altContext -> altContext.primaryReadContext().tumorQuality() >= config.filter().hardMinTumorQual() || tierSelector.isHotspot(
                altContext);
    }
}
