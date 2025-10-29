package com.hartwig.hmftools.sage.pipeline;

import static com.hartwig.hmftools.sage.SageCommon.SG_LOGGER;
import static com.hartwig.hmftools.sage.pipeline.ChromosomePartition.getPanelRegions;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.region.BaseRegion;
import com.hartwig.hmftools.common.region.ChrBaseRegion;
import com.hartwig.hmftools.sage.SageCallConfig;
import com.hartwig.hmftools.sage.candidate.Candidate;
import com.hartwig.hmftools.sage.candidate.Candidates;
import com.hartwig.hmftools.sage.candidate.ReadContextCandidate;
import com.hartwig.hmftools.sage.common.SamSlicerFactory;
import com.hartwig.hmftools.sage.common.SamSlicerInterface;
import com.hartwig.hmftools.common.variant.SimpleVariant;
import com.hartwig.hmftools.sage.evidence.CandidateEvidence;
import com.hartwig.hmftools.sage.common.RefSequence;

public class CandidateStage
{
    private final SageCallConfig mConfig;
    private final List<SimpleVariant> mHotspots;
    private final List<BaseRegion> mPanelRegions;
    private final CandidateEvidence mCandidateEvidence;
    private final List<BaseRegion> mHighConfidenceRegions;
    private final SamSlicerFactory mSamSlicerFactory;

    public CandidateStage(
            final SageCallConfig config, final List<SimpleVariant> hotspots,
            final List<BaseRegion> panelRegions, final List<BaseRegion> highConfidenceRegions,
            final SamSlicerFactory samSlicerFactory)
    {
        mConfig = config;
        mHotspots = hotspots;
        mPanelRegions = panelRegions;
        mHighConfidenceRegions = highConfidenceRegions;
        mSamSlicerFactory = samSlicerFactory;

        mCandidateEvidence = new CandidateEvidence(config.Common, hotspots, panelRegions);
    }

    public int totalReadsProcessed() { return mCandidateEvidence.totalReadsProcessed(); }

    public List<Candidate> findCandidates(final ChrBaseRegion region, final RefSequence refSequence)
    {
        final Candidates initialCandidates = new Candidates(mHotspots, mPanelRegions, mHighConfidenceRegions);

        final List<ChrBaseRegion> sliceRegions = !mConfig.PanelOnly ? Lists.newArrayList(region) : getPanelRegions(region, mPanelRegions);

        for(int i = 0; i < mConfig.TumorIds.size(); i++)
        {
            final String sample = mConfig.TumorIds.get(i);

            SamSlicerInterface samSlicer = mSamSlicerFactory.getSamSlicer(sample, sliceRegions, true);

            List<ReadContextCandidate> altCandidates = mCandidateEvidence.readBam(samSlicer, refSequence, region);

            if(mConfig.TumorIds.size() == 1)
                initialCandidates.addSingleSample(altCandidates);
            else
                initialCandidates.addOfMultipleSamples(altCandidates);
        }

        List<Candidate> candidates = initialCandidates.candidates(mConfig.Common.SpecificVariants);

        SG_LOGGER.trace("region({}) found {} candidates", region.toString(), candidates.size());

        return candidates;
    }
}
