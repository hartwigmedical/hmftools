package com.hartwig.hmftools.sage.evidence;

import static com.hartwig.hmftools.common.utils.sv.BaseRegion.positionWithin;
import static com.hartwig.hmftools.sage.evidence.ReadMatchType.NO_SUPPORT;
import static com.hartwig.hmftools.sage.evidence.ReadMatchType.SUPPORT;

import java.util.List;
import java.util.Map;

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

import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.reference.ReferenceSequenceFile;

public class ReadContextEvidence
{
    private final int mTypicalReadLength;
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

    public ReadContextEvidence(
            final SageConfig config, final RefGenomeInterface refGenome,
            final Map<String,QualityRecalibrationMap> qualityRecalibrationMap)
    {
        mSageConfig = config;
        mRefGenome = refGenome;
        mFactory = new ReadContextCounterFactory(config);
        mTypicalReadLength = config.typicalReadLength();
        mQualityRecalibrationMap = qualityRecalibrationMap;

        mRefSequence = null;
        mQualityCalculator = null;
        mReadCounters = null;
        mLastCandidateIndex = 0;
        mVariantPhaser = null;
    }

    public List<ReadContextCounter> collectEvidence(
            final List<Candidate> candidates, final String sample, final SamSlicerFactory samSlicerFactory, final VariantPhaser variantPhaser)
    {
        mReadCounters = mFactory.create(candidates);
        mLastCandidateIndex = 0;

        if(candidates.isEmpty())
            return mReadCounters;

        final Candidate firstCandidate = candidates.get(0);
        final Candidate lastCandidate = candidates.get(candidates.size() - 1);

        final ChrBaseRegion sliceRegion = new ChrBaseRegion(
                firstCandidate.chromosome(),
                Math.max(firstCandidate.position() - mTypicalReadLength, 1), lastCandidate.position() + mTypicalReadLength);

        mVariantPhaser = variantPhaser;

        if(mVariantPhaser != null)
            mVariantPhaser.initialise(sliceRegion, mSageConfig.LogLpsData);

        mRefSequence = new RefSequence(sliceRegion, mRefGenome);

        QualityRecalibrationMap qrMap = mQualityRecalibrationMap.get(sample);
        mQualityCalculator = new QualityCalculator(mSageConfig.Quality, qrMap, mRefSequence.IndexedBases);

        final SamSlicerInterface samSlicer = samSlicerFactory.getSamSlicer(sample, Lists.newArrayList(sliceRegion));
        samSlicer.slice(this::processReadRecord);

        return mReadCounters;
    }

    private void processReadRecord(final SAMRecord record)
    {
        // find any candidate potentially interested in this record
        int readStart = record.getAlignmentStart();
        int readEnd = record.getAlignmentEnd();

        if(record.getCigar().getFirstCigarElement().getOperator() == CigarOperator.S)
        {
            readStart -= record.getCigar().getFirstCigarElement().getLength();
        }

        if(record.getCigar().getLastCigarElement().getOperator() == CigarOperator.S)
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
            ReadMatchType matchType = readCounter.processRead(record, mSageConfig, mQualityCalculator, numberOfEvents);

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
