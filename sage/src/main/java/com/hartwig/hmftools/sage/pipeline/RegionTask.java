package com.hartwig.hmftools.sage.pipeline;

import static com.hartwig.hmftools.sage.SageCommon.SG_LOGGER;

import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.concurrent.Callable;
import java.util.function.Consumer;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.gene.TranscriptData;
import com.hartwig.hmftools.common.utils.PerformanceCounter;
import com.hartwig.hmftools.common.utils.sv.BaseRegion;
import com.hartwig.hmftools.common.utils.sv.ChrBaseRegion;
import com.hartwig.hmftools.common.variant.hotspot.VariantHotspot;
import com.hartwig.hmftools.sage.candidate.Candidate;
import com.hartwig.hmftools.sage.common.RefSequence;
import com.hartwig.hmftools.sage.common.SageVariant;
import com.hartwig.hmftools.sage.config.SageConfig;
import com.hartwig.hmftools.sage.config.VariantFilters;
import com.hartwig.hmftools.sage.coverage.Coverage;
import com.hartwig.hmftools.sage.evidence.ReadContextCounter;
import com.hartwig.hmftools.sage.evidence.ReadContextCounters;
import com.hartwig.hmftools.sage.evidence.VariantPhaser;
import com.hartwig.hmftools.sage.phase.PhaseSetCounter;
import com.hartwig.hmftools.sage.phase.VariantDeduper;
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
    private final VariantDeduper mVariantDeduper;

    private final List<SageVariant> mSageVariants;
    private final Set<Integer> mPassingPhaseSets;

    private final List<PerformanceCounter> mPerfCounters;

    public static final int PC_CANDIDATES = 0;
    public static final int PC_EVIDENCE = 1;
    public static final int PC_DEDUP = 2;

    public RegionTask(
            final int taskId, final ChrBaseRegion region, final SageConfig config, final ReferenceSequenceFile refGenome,
            final List<VariantHotspot> hotspots, final List<BaseRegion> panelRegions, final List<TranscriptData> transcripts,
            final List<BaseRegion> highConfidenceRegions, final Map<String, QualityRecalibrationMap> qualityRecalibrationMap,
            final PhaseSetCounter phaseSetCounter, final Coverage coverage)
    {
        mTaskId = taskId;
        mRegion = region;
        mConfig = config;
        mRefGenome = refGenome;

        mCandidateState = new CandidateStage(config, refGenome, hotspots, panelRegions, highConfidenceRegions, coverage);
        mEvidenceStage = new EvidenceStage(config, refGenome, qualityRecalibrationMap, phaseSetCounter);

        mVariantDeduper = new VariantDeduper(transcripts, phaseSetCounter, this::acceptDedupedVariant);

        mSageVariants = Lists.newArrayList();
        mPassingPhaseSets = Sets.newHashSet();

        mPerfCounters = Lists.newArrayList();
        mPerfCounters.add(new PerformanceCounter("Candidates"));
        mPerfCounters.add(new PerformanceCounter("Evidence"));
        // mPerfCounters.add(new PerformanceCounter("Dedup"));
    }

    public final List<SageVariant> getVariants() { return mSageVariants; }

    public final List<PerformanceCounter> getPerfCounters()
    {
        mPerfCounters.addAll(mEvidenceStage.getVariantPhaser().getPerfCounters());
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

        // final SageVariantFactory variantFactory = new SageVariantFactory(mConfig.Filter);
        VariantFilters filters = new VariantFilters(mConfig.Filter);

        // combine normal and tumor together to create variants, then apply soft filters
        Set<ReadContextCounter> passingTumorReadCounters = Sets.newHashSet();
        Set<ReadContextCounter> validTumorReadCounters = Sets.newHashSet(); // those not hard-filtered

        for(Candidate candidate : finalCandidates)
        {
            final List<ReadContextCounter> normalReadCounters = normalEvidence.getVariantReadCounters(candidate.variant());
            final List<ReadContextCounter> tumorReadCounters = tumorEvidence.getVariantReadCounters(candidate.variant());

            SageVariant sageVariant = new SageVariant(candidate, normalReadCounters, tumorReadCounters);
            mSageVariants.add(sageVariant);

            // apply filters
            if(filters.enabled())
                filters.applySoftFilters(sageVariant);

            if(sageVariant.isPassing())
                passingTumorReadCounters.add(tumorReadCounters.get(0));

            validTumorReadCounters.add(tumorReadCounters.get(0));
        }

        // phase variants now all evidence has been collected and filters applied
        VariantPhaser variantPhaser = mEvidenceStage.getVariantPhaser();

        variantPhaser.assignLocalPhaseSets(passingTumorReadCounters, validTumorReadCounters);

        // mPerfCounters.get(PC_DEDUP).start();

        SG_LOGGER.trace("phasing {} variants", mSageVariants.size());
        mSageVariants.forEach(mVariantDeduper);
        mVariantDeduper.flush();

        // mPerfCounters.get(PC_DEDUP).stop();

        SG_LOGGER.trace("{}: region({}) complete", mTaskId, mRegion);

        return (long)0;
    }

    private void acceptDedupedVariant(final SageVariant variant)
    {
        if(variant.isPassing() && variant.hasLocalPhaseSets())
            mPassingPhaseSets.addAll(variant.localPhaseSets());

        /*
        if(checkWriteVariant(variant, mVariantDeduper.passingPhaseSets()))
        {
            mFinalSageVariants.add(variant);
            // mWriteConsumer.accept(variant);
        }
        */
    }

    public void writeVariants(final Consumer<SageVariant> variantWriter)
    {
        List<SageVariant> finalVariants = mSageVariants.stream()
                .filter(x -> VariantFilters.checkFinalFilters(x, mPassingPhaseSets, mConfig)).collect(Collectors.toList());

        VariantPhaser.removeUninformativeLps(finalVariants, mPassingPhaseSets);

        finalVariants.forEach(x -> variantWriter.accept(x));
    }

}
