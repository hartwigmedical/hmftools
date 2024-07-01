package com.hartwig.hmftools.sage.pipeline;

import static java.lang.String.format;

import static com.hartwig.hmftools.sage.SageCommon.SG_LOGGER;
import static com.hartwig.hmftools.sage.common.RepeatInfo.setReferenceMaxRepeatInfo;

import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.gene.TranscriptData;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeInterface;
import com.hartwig.hmftools.common.region.BaseRegion;
import com.hartwig.hmftools.common.region.ChrBaseRegion;
import com.hartwig.hmftools.common.utils.PerformanceCounter;
import com.hartwig.hmftools.sage.SageCallConfig;
import com.hartwig.hmftools.sage.bqr.BqrRecordMap;
import com.hartwig.hmftools.sage.candidate.Candidate;
import com.hartwig.hmftools.sage.common.RefSequence;
import com.hartwig.hmftools.sage.common.SageVariant;
import com.hartwig.hmftools.sage.common.SamSlicerFactory;
import com.hartwig.hmftools.sage.common.SimpleVariant;
import com.hartwig.hmftools.sage.coverage.Coverage;
import com.hartwig.hmftools.sage.dedup.VariantDeduper;
import com.hartwig.hmftools.common.sage.FragmentLengthCounts;
import com.hartwig.hmftools.sage.evidence.FragmentLengths;
import com.hartwig.hmftools.sage.evidence.ReadContextCounter;
import com.hartwig.hmftools.sage.evidence.ReadContextCounters;
import com.hartwig.hmftools.sage.filter.VariantFilters;
import com.hartwig.hmftools.sage.phase.PhaseSetCounter;
import com.hartwig.hmftools.sage.phase.VariantPhaser;
import com.hartwig.hmftools.sage.quality.MsiJitterCalcs;
import com.hartwig.hmftools.sage.vis.VariantVis;

public class RegionTask
{
    private final ChrBaseRegion mRegion; // region to slice and analyse for this task
    private final int mTaskId;
    private final RegionResults mResults;
    private final FragmentLengths mFragmentLengths;

    private final SageCallConfig mConfig;
    private final RefGenomeInterface mRefGenome;

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
            final int taskId, final ChrBaseRegion region, final RegionResults results, final SageCallConfig config,
            final RefGenomeInterface refGenome, final List<SimpleVariant> hotspots, final List<BaseRegion> panelRegions,
            final List<TranscriptData> transcripts, final List<BaseRegion> highConfidenceRegions,
            final Map<String, BqrRecordMap> qualityRecalibrationMap, final MsiJitterCalcs msiJitterCalcs, final PhaseSetCounter phaseSetCounter,
            final Coverage coverage, final SamSlicerFactory samSlicerFactory, final FragmentLengths fragmentLengths)
    {
        mTaskId = taskId;
        mRegion = region;
        mResults = results;
        mConfig = config;
        mRefGenome = refGenome;
        mFragmentLengths = fragmentLengths;

        mCandidateState = new CandidateStage(config, hotspots, panelRegions, highConfidenceRegions, coverage, samSlicerFactory);

        mEvidenceStage = new EvidenceStage(
                config.Common, refGenome, qualityRecalibrationMap, msiJitterCalcs, phaseSetCounter, samSlicerFactory);

        mVariantDeduper = new VariantDeduper(transcripts, mRefGenome, mConfig.Common.getReadLength(), mConfig.Common.Filter);

        mSageVariants = Lists.newArrayList();
        mPassingPhaseSets = Sets.newHashSet();

        mPerfCounters = Lists.newArrayList();
        mPerfCounters.add(new PerformanceCounter("Candidates"));
        mPerfCounters.add(new PerformanceCounter("Evidence"));
        mPerfCounters.add(new PerformanceCounter("Variants"));
    }

    public final List<SageVariant> getVariants() { return mSageVariants; }

    public void run()
    {
        SG_LOGGER.trace("{}: region({}) finding candidates", mTaskId, mRegion);

        final RefSequence refSequence = new RefSequence(mRegion, mRefGenome);

        mPerfCounters.get(PC_CANDIDATES).start();
        List<Candidate> initialCandidates = mCandidateState.findCandidates(mRegion, refSequence);
        mPerfCounters.get(PC_CANDIDATES).stop();

        if(mConfig.Common.PerfWarnTime > 0 && mPerfCounters.get(PC_CANDIDATES).getLastTime() > mConfig.Common.PerfWarnTime)
        {
            SG_LOGGER.warn("region({}) candidate({}) reads({}) processing time({})",
                    mRegion, initialCandidates.size(), mCandidateState.totalReadsProcessed(),
                    String.format("%.3f", mPerfCounters.get(PC_CANDIDATES).getLastTime()));
        }

        if(initialCandidates.isEmpty())
        {
            SG_LOGGER.trace("{}: region({}) complete with no candidates", mTaskId, mRegion);
            return;
        }

        mResults.addCandidates(initialCandidates.size());

        SG_LOGGER.trace("{}: region({}) building evidence for {} candidates", mTaskId, mRegion, initialCandidates.size());

        mPerfCounters.get(PC_EVIDENCE).start();

        ReadContextCounters tumorEvidence = mEvidenceStage.findEvidence(
                mRegion, "tumor", mConfig.TumorIds, initialCandidates, true);

        List<Candidate> finalCandidates = tumorEvidence.filterCandidates();

        ReadContextCounters referenceEvidence = mEvidenceStage.findEvidence
                (mRegion, "reference", mConfig.Common.ReferenceIds, finalCandidates, false);

        mPerfCounters.get(PC_EVIDENCE).stop();

        VariantPhaser variantPhaser = mEvidenceStage.getVariantPhaser();

        if(mConfig.Common.PerfWarnTime > 0 && mPerfCounters.get(PC_EVIDENCE).getLastTime() > mConfig.Common.PerfWarnTime)
        {
            SG_LOGGER.warn("region({}) evidence candidates({}) phasing(g={} c={}) hardFilter({}) processing time({})",
                    mRegion, finalCandidates.size(),  variantPhaser.getPhasingGroupCount(), variantPhaser.getPhasedCollections().size(),
                    tumorEvidence.variantFilters().filterCountsStr(), String.format("%.3f", mPerfCounters.get(PC_EVIDENCE).getLastTime()));
        }

        variantPhaser.signalPhaseReadsEnd();

        if(!finalCandidates.isEmpty())
        {
            mPerfCounters.get(PC_VARIANTS).start();

            VariantFilters filters = new VariantFilters(mConfig.Common);

            // combine reference and tumor together to create variants, then apply soft filters
            Set<ReadContextCounter> passingTumorReadCounters = Sets.newHashSet();
            Set<ReadContextCounter> validTumorReadCounters = Sets.newHashSet(); // those not hard-filtered

            for(int candidateIndex = 0; candidateIndex < finalCandidates.size(); ++candidateIndex)
            {
                Candidate candidate = finalCandidates.get(candidateIndex);

                final List<ReadContextCounter> refCounters = !mConfig.Common.ReferenceIds.isEmpty() ?
                        referenceEvidence.getReadCounters(candidateIndex) : Lists.newArrayList();

                final List<ReadContextCounter> tumorReadCounters = tumorEvidence.getFilteredReadCounters(candidateIndex);

                SageVariant sageVariant = new SageVariant(candidate, refCounters, tumorReadCounters);
                setReferenceMaxRepeatInfo(sageVariant, refSequence);
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

            SG_LOGGER.trace("region({}) phasing {} variants", mRegion, mSageVariants.size());

            mVariantDeduper.processVariants(mSageVariants);

            mPerfCounters.get(PC_VARIANTS).stop();

            finaliseResults();
        }

        SG_LOGGER.trace("{}: region({}) complete", mTaskId, mRegion);
    }

    private void finaliseResults()
    {
        mSageVariants.stream().filter(x -> x.isPassing() && x.hasLocalPhaseSets()).forEach(x -> mPassingPhaseSets.addAll(x.localPhaseSets()));

        List<SageVariant> finalVariants = mSageVariants.stream()
                .filter(x -> VariantFilters.checkFinalFilters(x, mPassingPhaseSets, mConfig.Common, mConfig.PanelOnly))
                .collect(Collectors.toList());

        VariantPhaser.removeUninformativeLps(finalVariants, mPassingPhaseSets);

        mResults.addFinalVariants(mTaskId, finalVariants);

        if(mConfig.Common.Visualiser.Enabled)
        {
            mSageVariants.forEach(variant -> VariantVis.writeToHtmlFile(
                    variant, mConfig.TumorIds, mConfig.Common.ReferenceIds, mConfig.Common.Visualiser));
        }

        mResults.addTotalReads(mCandidateState.totalReadsProcessed());

        mPerfCounters.add(mEvidenceStage.getVariantPhaser().getPerfCounter());

        if(mConfig.Common.logPerfStats())
            mResults.addPerfCounters(mPerfCounters);

        mResults.addSynCounts(mEvidenceStage.getSyncCounts());
        mResults.addEvidenceStats(mEvidenceStage.getEvidenceStats());

        if(mConfig.Common.WriteFragmentLengths)
        {
            for(SageVariant variant : mSageVariants)
            {
                String variantInfo = format("%s:%d %s>%s", variant.chromosome(), variant.position(), variant.ref(), variant.alt());

                for(int s = 0; s < mConfig.TumorIds.size(); ++s)
                {
                    String sampleId = mConfig.TumorIds.get(s);
                    FragmentLengthCounts fragmentLengthData = variant.tumorReadCounters().get(s).fragmentLengths();
                    mFragmentLengths.writeVariantFragmentLength(variantInfo, sampleId, fragmentLengthData);
                }
            }
        }
    }
}
