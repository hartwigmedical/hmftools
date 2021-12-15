package com.hartwig.hmftools.sage.pipeline;

import static com.hartwig.hmftools.sage.SageCommon.SG_LOGGER;

import java.util.List;
import java.util.concurrent.CompletableFuture;

import com.hartwig.hmftools.common.utils.sv.ChrBaseRegion;
import com.hartwig.hmftools.common.variant.hotspot.VariantHotspot;
import com.hartwig.hmftools.sage.candidate.Candidate;
import com.hartwig.hmftools.sage.candidate.Candidates;
import com.hartwig.hmftools.sage.config.SageConfig;
import com.hartwig.hmftools.sage.context.AltContext;
import com.hartwig.hmftools.sage.coverage.Coverage;
import com.hartwig.hmftools.sage.evidence.CandidateEvidence;
import com.hartwig.hmftools.sage.ref.RefSequence;
import com.hartwig.hmftools.sage.sam.SamSlicerFactory;

import org.jetbrains.annotations.NotNull;

import htsjdk.samtools.reference.ReferenceSequenceFile;

public class CandidateStage
{
    private final SageConfig mConfig;
    private final List<VariantHotspot> mHotspots;
    private final List<ChrBaseRegion> mPanelRegions;
    private final CandidateEvidence mCandidateEvidence;
    private final List<ChrBaseRegion> mHighConfidenceRegions;

    public CandidateStage(final SageConfig config, final ReferenceSequenceFile refGenome,
            final List<VariantHotspot> hotspots, final List<ChrBaseRegion> panelRegions,
            final List<ChrBaseRegion> highConfidenceRegions, final Coverage coverage)
    {
        mConfig = config;

        final SamSlicerFactory samSlicerFactory = new SamSlicerFactory(config, panelRegions);

        mHotspots = hotspots;
        mPanelRegions = panelRegions;
        mHighConfidenceRegions = highConfidenceRegions;
        mCandidateEvidence = new CandidateEvidence(config, hotspots, panelRegions, samSlicerFactory, refGenome, coverage);
    }

    @NotNull
    public CompletableFuture<List<Candidate>> findCandidates(final ChrBaseRegion region, final CompletableFuture<RefSequence> refSequenceFuture)
    {
        return refSequenceFuture.thenCompose(refSequence ->
        {
            if(region.start() == 1)
            {
                SG_LOGGER.info("processing chromosome {}", region.Chromosome);
            }

            SG_LOGGER.trace("processing candidates in {}:{}", region.Chromosome, region.start());

            final Candidates initialCandidates = new Candidates(mHotspots, mPanelRegions, mHighConfidenceRegions);

            CompletableFuture<Void> done = CompletableFuture.completedFuture(null);

            for(int i = 0; i < mConfig.TumorIds.size(); i++)
            {
                final String sample = mConfig.TumorIds.get(i);
                final String sampleBam = mConfig.TumorBams.get(i);

                done = done.<List<AltContext>>thenApply(aVoid -> mCandidateEvidence.readBam(sample, sampleBam, refSequence, region))
                        .thenAccept(initialCandidates::add);
            }
            return done.thenApply(y -> initialCandidates.candidates());
        });
    }
}
