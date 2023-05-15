package com.hartwig.hmftools.sage.evidence;

import static java.lang.Math.abs;
import static java.lang.Math.max;
import static java.lang.Math.min;
import static java.lang.String.format;

import static com.hartwig.hmftools.common.utils.sv.BaseRegion.positionWithin;
import static com.hartwig.hmftools.sage.evidence.ReadMatchType.NO_SUPPORT;
import static com.hartwig.hmftools.sage.evidence.ReadMatchType.SUPPORT;

import static htsjdk.samtools.CigarOperator.S;

import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeInterface;
import com.hartwig.hmftools.common.utils.sv.ChrBaseRegion;
import com.hartwig.hmftools.sage.candidate.Candidate;
import com.hartwig.hmftools.sage.SageConfig;
import com.hartwig.hmftools.sage.common.SamSlicerFactory;
import com.hartwig.hmftools.sage.common.SamSlicerInterface;
import com.hartwig.hmftools.sage.phase.VariantPhaser;
import com.hartwig.hmftools.sage.quality.QualityCalculator;
import com.hartwig.hmftools.sage.quality.QualityRecalibrationMap;
import com.hartwig.hmftools.sage.common.RefSequence;
import com.hartwig.hmftools.sage.read.NumberEvents;

import htsjdk.samtools.SAMRecord;

public class ReadContextEvidence implements FragmentSyncReadHandler
{
    private final SageConfig mSageConfig;
    private final RefGenomeInterface mRefGenome;
    private final ReadContextCounterFactory mFactory;
    private final Map<String,QualityRecalibrationMap> mQualityRecalibrationMap;

    // state per slice region
    private RefSequence mRefSequence;
    private QualityCalculator mQualityCalculator;
    private List<ReadContextCounter> mReadCounters; // has one per candidate
    private int mLastCandidateIndex;

    private VariantPhaser mVariantPhaser;
    private final FragmentSync mFragmentSync;
    private String mCurrentSample;

    public ReadContextEvidence(
            final SageConfig config, final RefGenomeInterface refGenome, final Map<String,QualityRecalibrationMap> qualityRecalibrationMap)
    {
        mSageConfig = config;
        mRefGenome = refGenome;
        mFactory = new ReadContextCounterFactory(config);
        mQualityRecalibrationMap = qualityRecalibrationMap;
        mFragmentSync = new FragmentSync(this);

        mRefSequence = null;
        mQualityCalculator = null;
        mReadCounters = null;
        mLastCandidateIndex = 0;
        mVariantPhaser = null;
        mCurrentSample = null;
    }

    public List<ReadContextCounter> collectEvidence(
            final List<Candidate> candidates, final String sample, final SamSlicerFactory samSlicerFactory, final VariantPhaser variantPhaser)
    {
        mReadCounters = mFactory.create(candidates);
        mLastCandidateIndex = 0;

        if(candidates.isEmpty())
            return mReadCounters;


        List<Candidate> deleteCandidates = candidates.stream().filter(x -> x.variant().isDelete()).collect(Collectors.toList());

        for(ReadContextCounter readContextCounter : mReadCounters)
        {
            int maxCloseDel = deleteCandidates.stream()
                    .filter(x -> abs(x.position() - readContextCounter.position()) < 50)
                    .mapToInt(x -> x.variant().ref().length() - 1).max().orElse(0);

            if(maxCloseDel >= 5)
                readContextCounter.setMaxCandidateDeleteLength(maxCloseDel);
        }

        final Candidate firstCandidate = candidates.get(0);
        final Candidate lastCandidate = candidates.get(candidates.size() - 1);

        final ChrBaseRegion sliceRegion = new ChrBaseRegion(
                firstCandidate.chromosome(),
                max(firstCandidate.position() - mSageConfig.ExpectedReadLength, 1),
                lastCandidate.position() + mSageConfig.ExpectedReadLength);

        mVariantPhaser = variantPhaser;

        if(mVariantPhaser != null)
            mVariantPhaser.initialise(sliceRegion, mSageConfig.LogLpsData);

        mRefSequence = new RefSequence(sliceRegion, mRefGenome);

        QualityRecalibrationMap qrMap = mQualityRecalibrationMap.get(sample);
        mQualityCalculator = new QualityCalculator(mSageConfig.Quality, qrMap, mRefSequence.IndexedBases);

        mCurrentSample = sample;
        mFragmentSync.clear();

        final SamSlicerInterface samSlicer = samSlicerFactory.getSamSlicer(sample, Lists.newArrayList(sliceRegion), false);
        samSlicer.slice(this::processReadRecord);

        return mReadCounters;
    }

    public final int[] getSynCounts() { return mFragmentSync.getSynCounts(); }

    private void processReadRecord(final SAMRecord record)
    {
        processReadRecord(record, true);
    }

    @Override
    public void processReadRecord(final SAMRecord record, boolean checkSync)
    {
        if(checkSync && mSageConfig.SyncFragments)
        {
            if(mFragmentSync.handleOverlappingReads(record))
                return;
        }

        // find any candidate potentially interested in this record
        int readStart = record.getAlignmentStart();
        int readEnd = record.getAlignmentEnd();

        if(record.getCigar().getFirstCigarElement().getOperator() == S)
        {
            readStart -= record.getCigar().getFirstCigarElement().getLength();;
        }

        if(record.getCigar().getLastCigarElement().getOperator() == S)
        {
            readEnd += record.getCigar().getLastCigarElement().getLength();
        }

        // first look back from the last-used index
        List<ReadContextCounter> readCounters = Lists.newArrayList();

        // first check previous candidates starting with the current
        int nextIndex = mLastCandidateIndex + 1;
        int prevIndex = mLastCandidateIndex;

        while(prevIndex >= 0 && !mReadCounters.isEmpty())
        {
            ReadContextCounter readCounter = mReadCounters.get(prevIndex);

            if(positionWithin(readCounter.position(), readStart, readEnd))
            {
                mLastCandidateIndex = prevIndex;
                readCounters.add(0, readCounter);
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

            if(positionWithin(readCounter.position(), readStart, readEnd))
            {
                readCounters.add(readCounter);
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

        if(readCounters.isEmpty())
            return;

        List<ReadContextCounter> posPhasedCounters = mVariantPhaser != null ? Lists.newArrayList() : null;
        List<ReadContextCounter> negPhasedCounters = mVariantPhaser != null ? Lists.newArrayList() : null;

        int numberOfEvents = NumberEvents.calc(record, mRefSequence);

        for(ReadContextCounter readCounter : readCounters)
        {
            ReadMatchType matchType = readCounter.processRead(
                    record, mSageConfig, mQualityCalculator, numberOfEvents, mSageConfig.LogEvidenceReads ? mCurrentSample : null);

            if(mVariantPhaser != null)
            {
                if(matchType == SUPPORT)
                    posPhasedCounters.add(readCounter);
                else if(matchType == NO_SUPPORT)
                    negPhasedCounters.add(readCounter);
            }
        }

        if(mVariantPhaser != null)
            mVariantPhaser.registeredPhasedVariants(posPhasedCounters, negPhasedCounters);
    }
}
