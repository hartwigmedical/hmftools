package com.hartwig.hmftools.sage.pipeline;

import static java.util.concurrent.CompletableFuture.supplyAsync;

import static com.hartwig.hmftools.sage.SageCommon.SG_LOGGER;

import java.util.List;
import java.util.Map;
import java.util.concurrent.Callable;
import java.util.concurrent.CompletableFuture;
import java.util.concurrent.Executor;

import com.google.common.collect.Lists;
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
    private final ChrBaseRegion mRegion; // region to slice and analyse for this tasl
    private final int mTaskId;

    private final SageConfig mConfig;
    private final ReferenceSequenceFile mRefGenome;

    private final CandidateStage mCandidateState;
    private final EvidenceStage mEvidenceStage;

    private final List<SageVariant> mSageVariants;

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
    }

    // public final CandidateStage mCandidateState;
    // public final EvidenceStage mEvidenceStage;

    public final List<SageVariant> getVariants() { return mSageVariants; }

    @Override
    public Long call()
    {
        SG_LOGGER.trace("{}: region({}) finding candidates", mTaskId, mRegion);

        final RefSequence refSequence = new RefSequence(mRegion, mRefGenome);
        List<Candidate> initialCandidates = mCandidateState.findCandidates(mRegion, refSequence);

        SG_LOGGER.trace("{}: region({}) building evidence for {} candidates", mTaskId, mRegion, initialCandidates.size());

        ReadContextCounters tumorEvidence = mEvidenceStage.findEvidence(
                mRegion, "tumor", mConfig.TumorIds, mConfig.TumorBams, initialCandidates, true);

        // final List<Candidate> finalCandidates = filteredCandidates(tumorEvidence);
        List<Candidate> finalCandidates = tumorEvidence.candidates(mConfig.Filter.readContextFilter());

        ReadContextCounters normalEvidence = mEvidenceStage.findEvidence
                (mRegion, "normal", mConfig.ReferenceIds, mConfig.ReferenceBams, finalCandidates, false);

        final SageVariantFactory variantFactory = new SageVariantFactory(mConfig.Filter);

        // combine normal and tumor together and create variants
        for(Candidate candidate : finalCandidates)
        {
            final List<ReadContextCounter> normal = normalEvidence.readContextCounters(candidate.variant());
            final List<ReadContextCounter> tumor = tumorEvidence.readContextCounters(candidate.variant());

            SageVariant sageVariant = variantFactory.create(candidate, normal, tumor);
            mSageVariants.add(sageVariant);
        }

        SG_LOGGER.trace("{}: region({}) complete", mTaskId, mRegion);

        return (long)0;
    }

    /*
    public CompletableFuture<List<SageVariant>> findVariants(final ChrBaseRegion region)
    {
        final CompletableFuture<RefSequence> refSequenceFuture = supplyAsync(() -> new RefSequence(region, mRefGenome), mExecutor);

        final CompletableFuture<List<Candidate>> initialCandidates = mCandidateState.findCandidates(region, refSequenceFuture);

        final CompletableFuture<ReadContextCounters> tumorEvidence =
                mEvidenceStage.findEvidence(region, "tumor", mConfig.TumorIds, mConfig.TumorBams, initialCandidates, true);

        final CompletableFuture<List<Candidate>> finalCandidates = filteredCandidates(tumorEvidence);

        final CompletableFuture<ReadContextCounters> normalEvidence =
                mEvidenceStage.findEvidence(region, "normal", mConfig.ReferenceIds, mConfig.ReferenceBams, finalCandidates, false);

        return createSageVariants(finalCandidates, tumorEvidence, normalEvidence);
    }

    private CompletableFuture<List<SageVariant>> createSageVariants(
            final CompletableFuture<List<Candidate>> candidates, final CompletableFuture<ReadContextCounters> doneTumor,
            final CompletableFuture<ReadContextCounters> doneNormal)
    {
        return doneNormal.thenCombine(doneTumor, (normalCandidates, tumorCandidates) ->
        {
            // SG_LOGGER.trace("gathering evidence in {}:{}", region.Chromosome, region.start());

            final SageVariantFactory variantFactory = new SageVariantFactory(mConfig.Filter);

            // Combine normal and tumor together and create variants
            final List<SageVariant> result = Lists.newArrayList();
            for(Candidate candidate : candidates.join())
            {
                final List<ReadContextCounter> normal = normalCandidates.readContextCounters(candidate.variant());
                final List<ReadContextCounter> tumor = tumorCandidates.readContextCounters(candidate.variant());
                SageVariant sageVariant = variantFactory.create(candidate, normal, tumor);
                result.add(sageVariant);
            }

            return result;
        });
    }

    private CompletableFuture<List<Candidate>> filteredCandidates(final CompletableFuture<ReadContextCounters> tumorEvidence)
    {
        return tumorEvidence.thenApply(x -> x.candidates(mConfig.Filter.readContextFilter()));
    }
    */
}
