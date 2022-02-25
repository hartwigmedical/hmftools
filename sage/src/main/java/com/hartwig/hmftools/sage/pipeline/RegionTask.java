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
import com.hartwig.hmftools.sage.SageConfig;
import com.hartwig.hmftools.sage.filter.VariantFilters;
import com.hartwig.hmftools.sage.coverage.Coverage;
import com.hartwig.hmftools.sage.evidence.ReadContextCounter;
import com.hartwig.hmftools.sage.evidence.ReadContextCounters;
import com.hartwig.hmftools.sage.phase.VariantPhaser;
import com.hartwig.hmftools.sage.phase.PhaseSetCounter;
import com.hartwig.hmftools.sage.dedup.VariantDeduper;
import com.hartwig.hmftools.sage.quality.QualityRecalibrationMap;

import htsjdk.samtools.reference.ReferenceSequenceFile;

public class RegionTask implements Callable
{
    private final ChrBaseRegion mRegion; // region to slice and analyse for this task
    private final int mTaskId;
    private final RegionResults mResults;

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
    public static final int PC_VARIANTS = 2;

    public RegionTask(
            final int taskId, final ChrBaseRegion region, final RegionResults results, final SageConfig config, final ReferenceSequenceFile refGenome,
            final List<VariantHotspot> hotspots, final List<BaseRegion> panelRegions, final List<TranscriptData> transcripts,
            final List<BaseRegion> highConfidenceRegions, final Map<String, QualityRecalibrationMap> qualityRecalibrationMap,
            final PhaseSetCounter phaseSetCounter, final Coverage coverage)
    {
        mTaskId = taskId;
        mRegion = region;
        mResults = results;
        mConfig = config;
        mRefGenome = refGenome;

        mCandidateState = new CandidateStage(config, refGenome, hotspots, panelRegions, highConfidenceRegions, coverage);
        mEvidenceStage = new EvidenceStage(config, refGenome, qualityRecalibrationMap, phaseSetCounter);

        mVariantDeduper = new VariantDeduper(transcripts);

        mSageVariants = Lists.newArrayList();
        mPassingPhaseSets = Sets.newHashSet();

        mPerfCounters = Lists.newArrayList();
        mPerfCounters.add(new PerformanceCounter("Candidates"));
        mPerfCounters.add(new PerformanceCounter("Evidence"));
        mPerfCounters.add(new PerformanceCounter("Variants"));
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

        if(mConfig.PerfWarnTime > 0 && mPerfCounters.get(PC_CANDIDATES).getLastTime() > mConfig.PerfWarnTime)
        {
            SG_LOGGER.warn("region({}) candidate({}) reads({}) processing time({})",
                    mRegion, initialCandidates.size(), mCandidateState.totalReadsProcessed(),
                    String.format("%.3f", mPerfCounters.get(PC_CANDIDATES).getLastTime()));
        }

        SG_LOGGER.trace("{}: region({}) building evidence for {} candidates", mTaskId, mRegion, initialCandidates.size());

        mPerfCounters.get(PC_EVIDENCE).start();

        ReadContextCounters tumorEvidence = mEvidenceStage.findEvidence(
                mRegion, "tumor", mConfig.TumorIds, mConfig.TumorBams, initialCandidates, true);

        List<Candidate> finalCandidates = tumorEvidence.filterCandidates();

        ReadContextCounters normalEvidence = mEvidenceStage.findEvidence
                (mRegion, "normal", mConfig.ReferenceIds, mConfig.ReferenceBams, finalCandidates, false);

        mPerfCounters.get(PC_EVIDENCE).stop();

        VariantPhaser variantPhaser = mEvidenceStage.getVariantPhaser();

        if(mConfig.PerfWarnTime > 0 && mPerfCounters.get(PC_EVIDENCE).getLastTime() > mConfig.PerfWarnTime)
        {
            SG_LOGGER.warn("region({}) evidence candidates({}) phasing(g={} c={}) hardFilter({}) processing time({})",
                    mRegion, finalCandidates.size(),  variantPhaser.getPhasingGroupCount(), variantPhaser.getPhasedCollections().size(),
                    tumorEvidence.variantFilters().filterCountsStr(), String.format("%.3f", mPerfCounters.get(PC_EVIDENCE).getLastTime()));
        }

        variantPhaser.signalPhaseReadsEnd();

        mPerfCounters.get(PC_VARIANTS).start();

        VariantFilters filters = new VariantFilters(mConfig.Filter);

        // combine normal and tumor together to create variants, then apply soft filters
        Set<ReadContextCounter> passingTumorReadCounters = Sets.newHashSet();
        Set<ReadContextCounter> validTumorReadCounters = Sets.newHashSet(); // those not hard-filtered

        for(int candidateIndex = 0; candidateIndex < finalCandidates.size(); ++candidateIndex)
        {
            Candidate candidate = finalCandidates.get(candidateIndex);

            final List<ReadContextCounter> normalReadCounters = !mConfig.ReferenceIds.isEmpty() ?
                    normalEvidence.getReadCounters(candidateIndex) : Lists.newArrayList();

            final List<ReadContextCounter> tumorReadCounters = tumorEvidence.getFilteredReadCounters(candidateIndex);

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
        variantPhaser.assignLocalPhaseSets(passingTumorReadCounters, validTumorReadCounters);
        variantPhaser.clearAll();

        SG_LOGGER.trace("phasing {} variants", mSageVariants.size());

        mVariantDeduper.processVariants(mSageVariants);

        mPerfCounters.get(PC_VARIANTS).stop();

        finaliseResults();

        SG_LOGGER.trace("{}: region({}) complete", mTaskId, mRegion);

        return (long)0;
    }

    private void finaliseResults()
    {
        mSageVariants.stream().filter(x -> x.isPassing() && x.hasLocalPhaseSets()).forEach(x -> mPassingPhaseSets.addAll(x.localPhaseSets()));

        List<SageVariant> finalVariants = mSageVariants.stream()
                .filter(x -> VariantFilters.checkFinalFilters(x, mPassingPhaseSets, mConfig)).collect(Collectors.toList());

        VariantPhaser.removeUninformativeLps(finalVariants, mPassingPhaseSets);

        mResults.addFinalVariants(mTaskId, finalVariants);

        mResults.addTotalReads(mCandidateState.totalReadsProcessed());

        mPerfCounters.addAll(mEvidenceStage.getVariantPhaser().getPerfCounters());
        mResults.addPerfCounters(mPerfCounters);
    }

    public void writeVariants(final Consumer<SageVariant> variantWriter)
    {
        mSageVariants.stream().filter(x -> x.isPassing() && x.hasLocalPhaseSets()).forEach(x -> mPassingPhaseSets.addAll(x.localPhaseSets()));

        List<SageVariant> finalVariants = mSageVariants.stream()
                .filter(x -> VariantFilters.checkFinalFilters(x, mPassingPhaseSets, mConfig)).collect(Collectors.toList());

        VariantPhaser.removeUninformativeLps(finalVariants, mPassingPhaseSets);

        finalVariants.forEach(x -> variantWriter.accept(x));
    }

}
