package com.hartwig.hmftools.sage.pipeline;

import static com.hartwig.hmftools.sage.SageCommon.SG_LOGGER;
import static com.hartwig.hmftools.sage.pipeline.ChromosomePartition.getPanelRegions;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.utils.sv.BaseRegion;
import com.hartwig.hmftools.common.utils.sv.ChrBaseRegion;
import com.hartwig.hmftools.common.variant.hotspot.VariantHotspot;
import com.hartwig.hmftools.sage.candidate.Candidate;
import com.hartwig.hmftools.sage.candidate.Candidates;
import com.hartwig.hmftools.sage.SageConfig;
import com.hartwig.hmftools.sage.candidate.AltContext;
import com.hartwig.hmftools.sage.common.SamSlicerFactory;
import com.hartwig.hmftools.sage.common.SamSlicerInterface;
import com.hartwig.hmftools.sage.coverage.Coverage;
import com.hartwig.hmftools.sage.evidence.CandidateEvidence;
import com.hartwig.hmftools.sage.common.RefSequence;

import htsjdk.samtools.reference.ReferenceSequenceFile;

public class CandidateStage
{
    private final SageConfig mConfig;
    private final List<VariantHotspot> mHotspots;
    private final List<BaseRegion> mPanelRegions;
    private final CandidateEvidence mCandidateEvidence;
    private final List<BaseRegion> mHighConfidenceRegions;
    private final SamSlicerFactory mSamSlicerFactory;

    public CandidateStage(
            final SageConfig config, final ReferenceSequenceFile refGenome, final List<VariantHotspot> hotspots,
            final List<BaseRegion> panelRegions, final List<BaseRegion> highConfidenceRegions, final Coverage coverage,
            final SamSlicerFactory samSlicerFactory)
    {
        mConfig = config;
        mHotspots = hotspots;
        mPanelRegions = panelRegions;
        mHighConfidenceRegions = highConfidenceRegions;
        mSamSlicerFactory = samSlicerFactory;

        mCandidateEvidence = new CandidateEvidence(config, hotspots, panelRegions, coverage);
    }

    public int totalReadsProcessed() { return mCandidateEvidence.totalReadsProcessed(); }

    public List<Candidate> findCandidates(final ChrBaseRegion region, final RefSequence refSequence)
    {
        final Candidates initialCandidates = new Candidates(mHotspots, mPanelRegions, mHighConfidenceRegions);

        final List<ChrBaseRegion> sliceRegions = !mConfig.PanelOnly ? Lists.newArrayList(region) : getPanelRegions(region, mPanelRegions);

        for(int i = 0; i < mConfig.TumorIds.size(); i++)
        {
            final String sample = mConfig.TumorIds.get(i);

            SamSlicerInterface samSlicer = mSamSlicerFactory.getSamSlicer(sample, sliceRegions);

            List<AltContext> altContexts = mCandidateEvidence.readBam(sample, samSlicer, refSequence, region);

            if(mConfig.TumorIds.size() == 1)
                initialCandidates.addSingleSample(altContexts);
            else
                initialCandidates.addOfMultipleSamples(altContexts);
        }

        List<Candidate> candidates = initialCandidates.candidates(mConfig.SpecificPositions);

        SG_LOGGER.trace("region({}) found {} candidates", region.toString(), candidates.size());

        return candidates;
    }
}
