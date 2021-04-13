package com.hartwig.hmftools.sage.pipeline;

import java.util.List;
import java.util.concurrent.CompletableFuture;

import com.hartwig.hmftools.common.genome.region.GenomeRegion;
import com.hartwig.hmftools.common.variant.hotspot.VariantHotspot;
import com.hartwig.hmftools.sage.candidate.Candidate;
import com.hartwig.hmftools.sage.candidate.Candidates;
import com.hartwig.hmftools.sage.config.SageConfig;
import com.hartwig.hmftools.sage.coverage.Coverage;
import com.hartwig.hmftools.sage.evidence.CandidateEvidence;
import com.hartwig.hmftools.sage.ref.RefSequence;
import com.hartwig.hmftools.sage.sam.SamSlicerFactory;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

import htsjdk.samtools.reference.ReferenceSequenceFile;

public class CandidateStage {

    private static final Logger LOGGER = LogManager.getLogger(CandidateStage.class);

    private final SageConfig config;
    private final List<VariantHotspot> hotspots;
    private final List<GenomeRegion> panelRegions;
    private final CandidateEvidence candidateEvidence;
    private final List<GenomeRegion> highConfidenceRegions;

    public CandidateStage(@NotNull final SageConfig config, @NotNull final ReferenceSequenceFile refGenome,
            @NotNull final List<VariantHotspot> hotspots, @NotNull final List<GenomeRegion> panelRegions,
            @NotNull final List<GenomeRegion> highConfidenceRegions, final Coverage coverage) {
        this.config = config;
        final SamSlicerFactory samSlicerFactory = new SamSlicerFactory(config, panelRegions);
        this.hotspots = hotspots;
        this.panelRegions = panelRegions;
        this.highConfidenceRegions = highConfidenceRegions;
        this.candidateEvidence = new CandidateEvidence(config, hotspots, panelRegions, samSlicerFactory, refGenome, coverage);
    }

    @NotNull
    public CompletableFuture<List<Candidate>> candidates(@NotNull final GenomeRegion region,
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
}
