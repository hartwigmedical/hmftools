package com.hartwig.hmftools.sage.pipeline;

import static com.hartwig.hmftools.sage.SageCommon.SG_LOGGER;

import java.util.List;
import java.util.concurrent.CompletableFuture;

import com.hartwig.hmftools.common.utils.sv.BaseRegion;
import com.hartwig.hmftools.common.utils.sv.ChrBaseRegion;
import com.hartwig.hmftools.common.variant.hotspot.VariantHotspot;
import com.hartwig.hmftools.sage.candidate.Candidate;
import com.hartwig.hmftools.sage.candidate.Candidates;
import com.hartwig.hmftools.sage.config.SageConfig;
import com.hartwig.hmftools.sage.candidate.AltContext;
import com.hartwig.hmftools.sage.coverage.Coverage;
import com.hartwig.hmftools.sage.evidence.CandidateEvidence;
import com.hartwig.hmftools.sage.common.RefSequence;

import org.jetbrains.annotations.NotNull;

import htsjdk.samtools.reference.ReferenceSequenceFile;

public class CandidateStage
{
    private final SageConfig mConfig;
    private final List<VariantHotspot> mHotspots;
    private final List<BaseRegion> mPanelRegions;
    private final CandidateEvidence mCandidateEvidence;
    private final List<BaseRegion> mHighConfidenceRegions;

    public CandidateStage(final SageConfig config, final ReferenceSequenceFile refGenome,
            final List<VariantHotspot> hotspots, final List<BaseRegion> panelRegions,
            final List<BaseRegion> highConfidenceRegions, final Coverage coverage)
    {
        mConfig = config;
        mHotspots = hotspots;
        mPanelRegions = panelRegions;
        mHighConfidenceRegions = highConfidenceRegions;

        mCandidateEvidence = new CandidateEvidence(config, hotspots, panelRegions, refGenome, coverage);
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

            // SG_LOGGER.trace("region({}) finding candidates", region.toString());

            final Candidates initialCandidates = new Candidates(mHotspots, mPanelRegions, mHighConfidenceRegions);

            CompletableFuture<Void> done = CompletableFuture.completedFuture(null);

            for(int i = 0; i < mConfig.TumorIds.size(); i++)
            {
                final String sample = mConfig.TumorIds.get(i);
                final String sampleBam = mConfig.TumorBams.get(i);

                // SG_LOGGER.trace("region({}) finding candidates from tumor sample({})", region, sample);

                done = done.<List<AltContext>>thenApply(aVoid -> mCandidateEvidence.readBam(sample, sampleBam, refSequence, region))
                        .thenAccept(initialCandidates::add);
            }

            //return done.thenApply(y -> initialCandidates.candidates());
            return done.thenApply(y -> collectCandidates(region, initialCandidates));
        });
    }

    private List<Candidate> collectCandidates(final ChrBaseRegion region, final Candidates initialCandidates)
    {
        List<Candidate> candidates = initialCandidates.candidates(mConfig.SpecificPositions);
        SG_LOGGER.trace("region({}) found {} candidates", region.toString(), candidates.size());

        return candidates;
    }
}
