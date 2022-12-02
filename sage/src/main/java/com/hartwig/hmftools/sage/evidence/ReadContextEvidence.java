package com.hartwig.hmftools.sage.evidence;

import static java.lang.Math.max;

import static com.hartwig.hmftools.common.utils.sv.BaseRegion.positionWithin;
import static com.hartwig.hmftools.common.utils.sv.BaseRegion.positionsOverlap;
import static com.hartwig.hmftools.sage.evidence.ReadMatchType.NO_SUPPORT;
import static com.hartwig.hmftools.sage.evidence.ReadMatchType.SUPPORT;

import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
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

import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.cram.encoding.readfeatures.SoftClip;

public class ReadContextEvidence
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
    private int mMaxDeleteLength;

    private VariantPhaser mVariantPhaser;
    private String mCurrentSample;
    private final Map<String,SAMRecord> mCachedReads;

    public ReadContextEvidence(
            final SageConfig config, final RefGenomeInterface refGenome, final Map<String,QualityRecalibrationMap> qualityRecalibrationMap)
    {
        mSageConfig = config;
        mRefGenome = refGenome;
        mFactory = new ReadContextCounterFactory(config);
        mQualityRecalibrationMap = qualityRecalibrationMap;

        mRefSequence = null;
        mQualityCalculator = null;
        mReadCounters = null;
        mLastCandidateIndex = 0;
        mMaxDeleteLength = 0;
        mVariantPhaser = null;
        mCurrentSample = null;
        mCachedReads = Maps.newHashMap();
    }

    public List<ReadContextCounter> collectEvidence(
            final List<Candidate> candidates, final String sample, final SamSlicerFactory samSlicerFactory, final VariantPhaser variantPhaser)
    {
        mReadCounters = mFactory.create(candidates);
        mLastCandidateIndex = 0;

        if(candidates.isEmpty())
            return mReadCounters;

        mMaxDeleteLength = candidates.stream()
                .filter(x -> x.variant().isIndel())
                .mapToInt(x -> max(x.variant().ref().length() - x.variant().alt().length(), 0)).max().orElse(0);

        if(mMaxDeleteLength >= 5)
            mReadCounters.forEach(x -> x.setMaxCandidateDeleteLength(mMaxDeleteLength));

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
        mCachedReads.clear();

        final SamSlicerInterface samSlicer = samSlicerFactory.getSamSlicer(sample, Lists.newArrayList(sliceRegion), false);
        samSlicer.slice(this::processReadRecord);

        return mReadCounters;
    }

    private boolean handleOverlappingReads(final SAMRecord record)
    {
        final SAMRecord otherRecord = mCachedReads.get(record.getReadName());

        if(otherRecord != null)
        {
            final SAMRecord fragmentRecord = formFragmentRead(otherRecord, record);
            mCachedReads.remove(record.getReadName());

            if(fragmentRecord != null)
            {
                processReadRecord(fragmentRecord, false);
            }
            else
            {
                // process both reads if a consensus failed
                processReadRecord(otherRecord, false);
                processReadRecord(record, false);
            }
        }

        if(!record.getContig().equals(record.getMateReferenceName()))
            return false;

        if(!positionsOverlap(
                record.getAlignmentStart(), record.getAlignmentEnd(),
                record.getMateAlignmentStart(), record.getMateAlignmentStart() + record.getReadLength()))
        {
            return false;
        }

        // cache until the paired read arrives
        mCachedReads.put(record.getReadName(), record);
        return true;
    }

    private void processReadRecord(final SAMRecord record)
    {
        processReadRecord(record, true);
    }

    private void processReadRecord(final SAMRecord record, boolean checkSync)
    {
        if(checkSync && mSageConfig.SyncFragments)
        {
            if(handleOverlappingReads(record))
                return;
        }

        // find any candidate potentially interested in this record
        int readStart = record.getAlignmentStart();
        int readEnd = record.getAlignmentEnd();

        if(record.getCigar().getFirstCigarElement().getOperator() == CigarOperator.S)
        {
            readStart -= record.getCigar().getFirstCigarElement().getLength();;
            readStart -= mMaxDeleteLength; // account for deleted bases being the cause of the soft-clipping
        }

        if(record.getCigar().getLastCigarElement().getOperator() == CigarOperator.S)
        {
            readEnd += record.getCigar().getLastCigarElement().getLength();
            readEnd += mMaxDeleteLength;
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

    private SAMRecord formFragmentRead(final SAMRecord first, final SAMRecord second)
    {
        // take the highest base qual base for any overlapping bases
        // widen the read to cover both reads
        // how to handle preserve INDELs?

        /*

        int firstPosStart = first.getAlignmentStart();
        int firstPosEnd = first.getAlignmentEnd();
        Cigar firstCigar = first.getCigar();
        final byte[] firstBaseQualities = first.getBaseQualities();
        final byte[] firstBases = first.getReadBases();

        int[] firstScLengths = new int[] { firstCigar.}

        int secondPosStart = second.getAlignmentStart();
        int secondPosEnd = second.getAlignmentEnd();
        Cigar secondCigar = second.getCigar();
        final byte[] secondBaseQualities = second.getBaseQualities();
        final byte[] secondBases = second.getReadBases();

        Cigar combinedCigar = new Cigar();
        final byte[] combinedBaseQualities = new byte[];
        final byte[] combinedBases = new byte[];
        int combinedPosStart;
        int combinedPosEnd;
        */

        return null;
    }
}
