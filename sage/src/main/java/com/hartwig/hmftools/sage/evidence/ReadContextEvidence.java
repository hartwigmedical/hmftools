package com.hartwig.hmftools.sage.evidence;

import static java.lang.Math.abs;
import static java.lang.Math.max;
import static java.lang.Math.min;

import static com.hartwig.hmftools.common.bam.SamRecordUtils.readToString;
import static com.hartwig.hmftools.common.region.BaseRegion.positionWithin;
import static com.hartwig.hmftools.common.region.BaseRegion.positionsOverlap;
import static com.hartwig.hmftools.common.region.BaseRegion.positionsWithin;
import static com.hartwig.hmftools.sage.SageCommon.SG_LOGGER;
import static com.hartwig.hmftools.sage.evidence.ReadMatchType.REF_SUPPORT;
import static com.hartwig.hmftools.sage.evidence.ReadMatchType.ALT_SUPPORT;

import static htsjdk.samtools.CigarOperator.N;

import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import javax.annotation.Nullable;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.bam.CigarUtils;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeInterface;
import com.hartwig.hmftools.common.region.ChrBaseRegion;
import com.hartwig.hmftools.sage.candidate.Candidate;
import com.hartwig.hmftools.sage.SageConfig;
import com.hartwig.hmftools.sage.common.SamSlicerFactory;
import com.hartwig.hmftools.sage.common.SamSlicerInterface;
import com.hartwig.hmftools.sage.phase.VariantPhaser;
import com.hartwig.hmftools.sage.quality.MsiJitterCalcs;
import com.hartwig.hmftools.sage.quality.QualityCalculator;
import com.hartwig.hmftools.sage.bqr.BqrRecordMap;
import com.hartwig.hmftools.sage.common.RefSequence;
import com.hartwig.hmftools.sage.common.NumberEvents;
import com.hartwig.hmftools.sage.sync.FragmentData;
import com.hartwig.hmftools.sage.sync.FragmentSync;
import com.hartwig.hmftools.sage.sync.FragmentSyncReadHandler;

import htsjdk.samtools.SAMRecord;

public class ReadContextEvidence implements FragmentSyncReadHandler
{
    private final SageConfig mConfig;
    private final RefGenomeInterface mRefGenome;
    private final ReadContextCounterFactory mFactory;
    private final Map<String,BqrRecordMap> mQualityRecalibrationMap;
    private final MsiJitterCalcs mMsiJitterCalcs;

    // state per slice region
    private RefSequence mRefSequence;
    private List<ReadContextCounter> mReadCounters; // has one per candidate
    private int mLastCandidateIndex;

    private VariantPhaser mVariantPhaser;
    private final FragmentSync mFragmentSync;
    private List<ReadContextCounter> mSelectedReadCounters; // persisted array to avoid re-creation on every read
    private final EvidenceStats mStats;

    public ReadContextEvidence(
            final SageConfig config, final RefGenomeInterface refGenome, final Map<String, BqrRecordMap> qualityRecalibrationMap,
            final MsiJitterCalcs msiJitterCalcs)
    {
        mConfig = config;
        mRefGenome = refGenome;
        mFactory = new ReadContextCounterFactory(config);
        mQualityRecalibrationMap = qualityRecalibrationMap;
        mMsiJitterCalcs = msiJitterCalcs;
        mFragmentSync = new FragmentSync(this, refGenome);

        mRefSequence = null;
        mReadCounters = null;
        mLastCandidateIndex = 0;
        mVariantPhaser = null;
        mSelectedReadCounters = null;
        mStats = new EvidenceStats();
    }

    private static final int SLICE_SOFT_CLIP_BUFFER = 30;

    public List<ReadContextCounter> collectEvidence(
            final List<Candidate> candidates, final String sample, final SamSlicerFactory samSlicerFactory, final VariantPhaser variantPhaser)
    {
        if(candidates.isEmpty())
            return Collections.emptyList();

        List<ChrBaseRegion> sliceRegions = buildCandidateRegions(candidates);

        // add a buffer around each slice to support variants in soft-clip regions
        for(ChrBaseRegion sliceRegion : sliceRegions)
        {
            sliceRegion.setStart(max(sliceRegion.start() - SLICE_SOFT_CLIP_BUFFER, 1));
            sliceRegion.setEnd(sliceRegion.end() + SLICE_SOFT_CLIP_BUFFER);
        }

        ++mStats.PartitionCount;
        mStats.SliceCount += sliceRegions.size();
        mStats.SliceLength += sliceRegions.stream().mapToInt(x -> x.baseLength()).sum();

        int sliceRegionStart = sliceRegions.get(0).start();
        int sliceRegionEnd = sliceRegions.get(sliceRegions.size() - 1).end();
        ChrBaseRegion regionBounds = new ChrBaseRegion(sliceRegions.get(0).chromosome(), sliceRegionStart, sliceRegionEnd);

        mVariantPhaser = variantPhaser;

        if(mVariantPhaser != null)
            mVariantPhaser.initialise(regionBounds, mConfig.LogLpsData);

        mRefSequence = new RefSequence(regionBounds, mRefGenome);

        BqrRecordMap qrMap = mQualityRecalibrationMap.get(sample);
        QualityCalculator qualityCalculator = new QualityCalculator(mConfig, qrMap, mRefSequence, mRefGenome, mMsiJitterCalcs);

        mReadCounters = mFactory.create(candidates, mConfig, qualityCalculator, sample);
        mLastCandidateIndex = 0;

        mSelectedReadCounters = Lists.newArrayListWithCapacity(mReadCounters.size());

        List<Candidate> deleteCandidates = candidates.stream().filter(x -> x.variant().isDelete()).collect(Collectors.toList());

        for(ReadContextCounter readContextCounter : mReadCounters)
        {
            int maxCloseDel = deleteCandidates.stream()
                    .filter(x -> abs(x.position() - readContextCounter.position()) < 50)
                    .mapToInt(x -> x.variant().ref().length() - 1).max().orElse(0);

            if(maxCloseDel >= 5)
                readContextCounter.setMaxCandidateDeleteLength(maxCloseDel);
        }

        final SamSlicerInterface samSlicer = samSlicerFactory.getSamSlicer(sample, sliceRegions, false);
        samSlicer.slice(this::processReadRecord);

        mFragmentSync.emptyCachedReads();

        if(mConfig.PerfWarnTime > 0)
        {
            SG_LOGGER.trace("region({}) evidence stats: {}", regionBounds, mStats);
        }

        if(mConfig.Quality.MapQualityRatioFactor > 0)
        {
            mReadCounters.forEach(x -> x.applyMapQualityRatio());
        }

        mReadCounters.forEach(x -> x.jitter().setJitterQualFilterState(mMsiJitterCalcs, x));

        return mReadCounters;
    }

    private List<ChrBaseRegion> buildCandidateRegions(final List<Candidate> candidates)
    {
        // make up to X regions based on distance between them
        final Candidate firstCandidate = candidates.get(0);
        final Candidate lastCandidate = candidates.get(candidates.size() - 1);

        List<ChrBaseRegion> sliceRegions = Lists.newArrayList();

        int minGap = max(mConfig.getReadLength() * 2, 100);

        int sliceLength = lastCandidate.position() - firstCandidate.position();
        int averageGap = sliceLength / candidates.size();

        if(averageGap <= minGap || mConfig.MaxPartitionSlices == 1)
        {
            sliceRegions.add(new ChrBaseRegion(firstCandidate.chromosome(), firstCandidate.position(), lastCandidate.position()));
            return sliceRegions;
        }

        List<Integer> gapDistances = Lists.newArrayListWithCapacity(candidates.size() - 1);

        for(int i = 0; i < candidates.size() - 1; ++i)
        {
            int gap = candidates.get(i + 1).position() - candidates.get(i).position();
            gapDistances.add(gap);
        }

        Collections.sort(gapDistances, Collections.reverseOrder());

        int nth = min(mConfig.MaxPartitionSlices, gapDistances.size());
        int nthGap = gapDistances.get(nth - 1);

        /*
        int gapCount = gapDistances.size();
        int medianGap = gapDistances.get(gapCount / 2);

        SG_LOGGER.debug("region({}:{}-{} len={}) candidates({}) gap(n={} nth={}, max={} avg={} median={})",
                firstCandidate.chromosome(), firstCandidate.position(), lastCandidate.position(), sliceLength,
                candidates.size(), nth, nthGap, gapDistances.get(0), averageGap, medianGap);
        */

        ChrBaseRegion currentRegion = new ChrBaseRegion(firstCandidate.chromosome(), firstCandidate.position(), firstCandidate.position());

        sliceRegions.add(currentRegion);

        for(int i = 1; i < candidates.size(); ++i)
        {
            Candidate candidate = candidates.get(i);
            int gap = candidate.position() - currentRegion.end();

            if(gap <= max(minGap, nthGap))
            {
                currentRegion.setEnd(candidate.position());
            }
            else
            {
                currentRegion = new ChrBaseRegion(candidate.chromosome(), candidate.position(), candidate.position());
                sliceRegions.add(currentRegion);
            }
        }

        return sliceRegions;
    }

    public final int[] getSynCounts() { return mFragmentSync.getSynCounts(); }
    public EvidenceStats evidenceStats() { return mStats; }

    private void processReadRecord(final SAMRecord record)
    {
        processReadRecord(record, true);
    }

    @Override
    public void processReadRecord(final SAMRecord record, boolean checkSync) { processReadRecord(record, checkSync, null); }

    @Override
    public void processReadRecord(final SAMRecord record, boolean checkSync, @Nullable final FragmentData fragmentData)
    {
        ++mStats.ReadCount;

        if(mConfig.SyncFragments && checkSync)
        {
            if(mFragmentSync.handleOverlappingReads(record))
                return;
        }

        mSelectedReadCounters.clear();

        // find any candidate potentially interested in this record
        int readStart = record.getAlignmentStart() - CigarUtils.leftSoftClipLength(record);
        int readEnd = record.getAlignmentEnd() + CigarUtils.rightSoftClipLength(record);

        // first check previous candidates starting with the current
        int nextIndex = mLastCandidateIndex + 1;
        int prevIndex = mLastCandidateIndex;

        while(prevIndex >= 0 && !mReadCounters.isEmpty())
        {
            ReadContextCounter readCounter = mReadCounters.get(prevIndex);

            if(considerRead(readCounter, readStart, readEnd))
            {
                mLastCandidateIndex = prevIndex;
                mSelectedReadCounters.add(0, readCounter);
                readStart -= readCounter.maxCandidateDeleteLength();
                readEnd += readCounter.maxCandidateDeleteLength();
            }
            else if(readCounter.position() < readStart)
            {
                break;
            }

            --prevIndex;
        }

        // now check ahead from the current index
        while(nextIndex < mReadCounters.size())
        {
            ReadContextCounter readCounter = mReadCounters.get(nextIndex);

            if(considerRead(readCounter, readStart, readEnd))
            {
                mSelectedReadCounters.add(readCounter);
            }
            else if(readCounter.position() < readStart)
            {
                mLastCandidateIndex = nextIndex;
            }
            else if(readCounter.position() > readEnd)
            {
                break;
            }

            ++nextIndex;
        }

        if(mSelectedReadCounters.isEmpty())
        {
            ++mStats.NoVariantReadCount;
            return;
        }

        List<ReadContextCounter> posPhasedCounters = mVariantPhaser != null ? Lists.newArrayList() : null;
        List<ReadContextCounter> negPhasedCounters = mVariantPhaser != null ? Lists.newArrayList() : null;

        int numberOfEvents = NumberEvents.calc(record, mRefSequence);

        for(ReadContextCounter readCounter : mSelectedReadCounters)
        {
            ReadMatchType matchType = null;

            try
            {
                matchType = readCounter.processRead(record, numberOfEvents, fragmentData);
            }
            catch(Exception e)
            {
                SG_LOGGER.error("var({}) read({}) error: {}", readCounter.readContext(), readToString(record), e.toString());
                e.printStackTrace();
                System.exit(1);
            }

            ++mStats.SupportCounts[matchType.ordinal()];

            if(mVariantPhaser != null)
            {
                if(matchType == ALT_SUPPORT)
                    posPhasedCounters.add(readCounter);
                else if(matchType == REF_SUPPORT)
                    negPhasedCounters.add(readCounter);
            }
        }

        if(mVariantPhaser != null)
        {
            mVariantPhaser.registeredPhasedVariants(posPhasedCounters, negPhasedCounters);
        }
    }

    private boolean considerRead(final ReadContextCounter readCounter, final int unclippedStart, final int unclippedEnd)
    {
        if(readCounter.exceedsMaxCoverage())
            return false;

        if(readCounter.variant().isDelete())
        {
            return positionsOverlap(readCounter.position(), readCounter.maxPositionVsReadStart(), unclippedStart, unclippedEnd);
        }

        return positionWithin(readCounter.position(), unclippedStart, unclippedEnd);
    }
}
