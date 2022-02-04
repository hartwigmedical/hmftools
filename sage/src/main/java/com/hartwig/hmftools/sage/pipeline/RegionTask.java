package com.hartwig.hmftools.sage.pipeline;

import static com.hartwig.hmftools.sage.SageCommon.SG_LOGGER;

import java.util.List;
import java.util.Map;
import java.util.concurrent.Callable;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.utils.PerformanceCounter;
import com.hartwig.hmftools.common.utils.sv.BaseRegion;
import com.hartwig.hmftools.common.utils.sv.ChrBaseRegion;
import com.hartwig.hmftools.common.variant.hotspot.VariantHotspot;
import com.hartwig.hmftools.sage.candidate.Candidate;
import com.hartwig.hmftools.sage.common.RefSequence;
import com.hartwig.hmftools.sage.common.SageVariant;
import com.hartwig.hmftools.sage.common.SageVariantFactory;
import com.hartwig.hmftools.sage.config.SageConfig;
import com.hartwig.hmftools.sage.coverage.Coverage;
import com.hartwig.hmftools.sage.evidence.ReadContextCounter;
import com.hartwig.hmftools.sage.evidence.ReadContextCounters;
import com.hartwig.hmftools.sage.phase.PhaseSetCounter;
import com.hartwig.hmftools.sage.quality.QualityRecalibrationMap;

import htsjdk.samtools.reference.ReferenceSequenceFile;

public class RegionTask implements Callable
{
    private final ChrBaseRegion mRegion; // region to slice and analyse for this task
    private final int mTaskId;

    private final SageConfig mConfig;
    private final ReferenceSequenceFile mRefGenome;

    private final CandidateStage mCandidateState;
    private final EvidenceStage mEvidenceStage;

    private final List<SageVariant> mSageVariants;

    private final List<PerformanceCounter> mPerfCounters;

    public static final int PC_CANDIDATES = 0;
    public static final int PC_EVIDENCE = 1;
    public static final int PC_DEDUP = 2;

    public RegionTask(
            final int taskId, final ChrBaseRegion region, final SageConfig config, final ReferenceSequenceFile refGenome,
            final List<VariantHotspot> hotspots, final List<BaseRegion> panelRegions,
            final List<BaseRegion> highConfidenceRegions, final Map<String, QualityRecalibrationMap> qualityRecalibrationMap,
            final PhaseSetCounter phaseSetCounter, final Coverage coverage)
    {
        mTaskId = taskId;
        mRegion = region;
        mConfig = config;
        mRefGenome = refGenome;

        mCandidateState = new CandidateStage(config, refGenome, hotspots, panelRegions, highConfidenceRegions, coverage);
        mEvidenceStage = new EvidenceStage(config, refGenome, qualityRecalibrationMap, phaseSetCounter);

        mSageVariants = Lists.newArrayList();

        mPerfCounters = Lists.newArrayList();
        mPerfCounters.add(new PerformanceCounter("Candidates"));
        mPerfCounters.add(new PerformanceCounter("Evidence"));
    }

    public final List<SageVariant> getVariants() { return mSageVariants; }

    public final List<PerformanceCounter> getPerfCounters()
    {
        mPerfCounters.addAll(mEvidenceStage.getPerfCounters());
        return mPerfCounters;
    }

    @Override
    public Long call()
    {
        SG_LOGGER.trace("{}: region({}) finding candidates", mTaskId, mRegion);

        final RefSequence refSequence = new RefSequence(mRegion, mRefGenome);

        mPerfCounters.get(PC_CANDIDATES).start();
        List<Candidate> initialCandidates = mCandidateState.findCandidates(mRegion, refSequence);
        mPerfCounters.get(PC_CANDIDATES).stop();

        SG_LOGGER.trace("{}: region({}) building evidence for {} candidates", mTaskId, mRegion, initialCandidates.size());

        mPerfCounters.get(PC_EVIDENCE).start();

        ReadContextCounters tumorEvidence = mEvidenceStage.findEvidence(
                mRegion, "tumor", mConfig.TumorIds, mConfig.TumorBams, initialCandidates, true);

        List<Candidate> finalCandidates = tumorEvidence.filterCandidates(mConfig.Filter);

        ReadContextCounters normalEvidence = mEvidenceStage.findEvidence
                (mRegion, "normal", mConfig.ReferenceIds, mConfig.ReferenceBams, finalCandidates, false);

        mPerfCounters.get(PC_EVIDENCE).stop();

        final SageVariantFactory variantFactory = new SageVariantFactory(mConfig.Filter);

        // combine normal and tumor together and create variants
        for(Candidate candidate : finalCandidates)
        {
            final List<ReadContextCounter> normal = normalEvidence.getVariantReadCounters(candidate.variant());
            final List<ReadContextCounter> tumor = tumorEvidence.getVariantReadCounters(candidate.variant());

            SageVariant sageVariant = variantFactory.create(candidate, normal, tumor);
            mSageVariants.add(sageVariant);
        }

        SG_LOGGER.trace("{}: region({}) complete", mTaskId, mRegion);

        return (long)0;
    }
}
