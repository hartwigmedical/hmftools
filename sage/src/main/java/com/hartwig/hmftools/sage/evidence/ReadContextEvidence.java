package com.hartwig.hmftools.sage.evidence;

import static com.hartwig.hmftools.common.utils.sv.BaseRegion.positionWithin;

import java.io.File;
import java.util.List;
import java.util.Map;
import java.util.Set;

import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.utils.sv.ChrBaseRegion;
import com.hartwig.hmftools.sage.candidate.Candidate;
import com.hartwig.hmftools.sage.config.SageConfig;
import com.hartwig.hmftools.sage.phase.PhaseSetCounter;
import com.hartwig.hmftools.sage.quality.QualityCalculator;
import com.hartwig.hmftools.sage.quality.QualityRecalibrationMap;
import com.hartwig.hmftools.sage.common.RefSequence;
import com.hartwig.hmftools.sage.common.SamSlicer;
import com.hartwig.hmftools.sage.read.NumberEvents;

import org.apache.commons.compress.utils.Lists;

import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.cram.ref.ReferenceSource;
import htsjdk.samtools.reference.ReferenceSequenceFile;

public class ReadContextEvidence
{
    private final int mTypicalReadLength;
    private final SageConfig mSageConfig;
    private final ReferenceSequenceFile mRefGenome;
    private final ReadContextCounterFactory mFactory;
    private final Map<String,QualityRecalibrationMap> mQualityRecalibrationMap;
    private final PhaseSetCounter mPhaseSetCounter;

    // initiated per slice region
    private RefSequence mRefSequence;
    private QualityCalculator mQualityCalculator;
    private List<ReadContextCounter> mReadCounters;
    private int mLastCandidateIndex;

    private final List<PhasedReadCounters> mPhasedReadCounters;

    public ReadContextEvidence(
            final SageConfig config, final ReferenceSequenceFile refGenome,
            final Map<String,QualityRecalibrationMap> qualityRecalibrationMap, final PhaseSetCounter phaseSetCounter)
    {
        mSageConfig = config;
        mRefGenome = refGenome;
        mFactory = new ReadContextCounterFactory(config);
        mTypicalReadLength = config.typicalReadLength();
        mQualityRecalibrationMap = qualityRecalibrationMap;
        mPhaseSetCounter = phaseSetCounter;

        mRefSequence = null;
        mQualityCalculator = null;
        mReadCounters = null;
        mLastCandidateIndex = 0;
        mPhasedReadCounters = Lists.newArrayList();
    }

    public List<ReadContextCounter> collectEvidence(final List<Candidate> candidates, final String sample, final String bam)
    {
        mReadCounters = mFactory.create(sample, candidates);

        if(candidates.isEmpty())
            return mReadCounters;

        final Candidate firstCandidate = candidates.get(0);
        final Candidate lastCandidate = candidates.get(candidates.size() - 1);

        final ChrBaseRegion bounds = new ChrBaseRegion(
                firstCandidate.chromosome(),
                Math.max(firstCandidate.position() - mTypicalReadLength, 1), lastCandidate.position() + mTypicalReadLength);

        mRefSequence = new RefSequence(bounds, mRefGenome);

        QualityRecalibrationMap qrMap = mQualityRecalibrationMap.get(sample);
        mQualityCalculator = new QualityCalculator(mSageConfig.Quality, qrMap, mRefSequence.IndexedBases);

        final SamReader tumorReader = SamReaderFactory.makeDefault().validationStringency(mSageConfig.Stringency)
                .referenceSource(new ReferenceSource(mRefGenome)).open(new File(bam));

        final SamSlicer slicer = new SamSlicer(0, bounds);

        mLastCandidateIndex = 0;
        slicer.slice(tumorReader, this::processReadRecord);

        /*

        try
        {
            final SamRecordSelector<ReadContextCounter> consumerSelector = new SamRecordSelector<>(counters);

            final RefSequence refSequence = new RefSequence(bounds, mRefGenome);

            QualityRecalibrationMap qrMap = mQualityRecalibrationMap.get(sample);

            mQualityCalculator = new QualityCalculator(mSageConfig.Quality, qrMap, refSequence.IndexedBases);

            final SamReader tumorReader = SamReaderFactory.makeDefault().validationStringency(mSageConfig.Stringency)
                    .referenceSource(new ReferenceSource(mRefGenome)).open(new File(bam));

            final SamSlicer slicer = new SamSlicer(0, bounds);

            slicer.slice(tumorReader, samRecord ->
            {
                int numberOfEvents = NumberEvents.numberOfEvents(samRecord, refSequence);
                consumerSelector.select(samRecord, x -> x.processRead(samRecord, mSageConfig, qualityCalculator, numberOfEvents));

            });
        }
        catch(IOException e)
        {
            throw new CompletionException(e);
        }
         */

        // assign local phase sets
        for(PhasedReadCounters phasedReadCounters : mPhasedReadCounters)
        {
            int nextLps = mPhaseSetCounter.getNext();
            phasedReadCounters.ReadCounters.forEach(x -> x.setLocalPhaseSet(nextLps));
        }

        mPhasedReadCounters.clear();

        return mReadCounters;
    }

    private void processReadRecord(final SAMRecord record)
    {
        int numberOfEvents = NumberEvents.numberOfEvents(record, mRefSequence);

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

        while(prevIndex >= 0)
        {
            ReadContextCounter readCounter = mReadCounters.get(prevIndex);

            // remove any counters at their limit for coverage
            while(readCounter.exceedsMaxCoverage())
            {
                mReadCounters.remove(prevIndex);
            }

            if(positionWithin(readCounter.position(), readStart, readEnd))
            {
                mLastCandidateIndex = prevIndex;
                readCounters.add(readCounter);
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

        Set<ReadContextCounter> phasedCounters = Sets.newHashSet();

        for(ReadContextCounter readCounter : readCounters)
        {
            if(readCounter.processRead(record, mSageConfig, mQualityCalculator, numberOfEvents))
                phasedCounters.add(readCounter);
        }

        if(phasedCounters.size() >= 2)
        {
            boolean matched = false;

            for(PhasedReadCounters phasedReadCounters : mPhasedReadCounters)
            {
                if(phasedReadCounters.processPhasedCounters(phasedCounters))
                {
                    matched = true;
                    break;
                }
            }

            if(!matched)
            {
                mPhasedReadCounters.add(new PhasedReadCounters(phasedCounters));
            }
        }
    }

    private class PhasedReadCounters
    {
        public final Set<ReadContextCounter> ReadCounters;
        public int ReadCount;

        public PhasedReadCounters(final Set<ReadContextCounter> phasedCounters)
        {
            ReadCounters = phasedCounters;
            ++ReadCount;
        }

        public boolean processPhasedCounters(final Set<ReadContextCounter> phasedCounters)
        {
            if(ReadCounters.stream().anyMatch(x -> phasedCounters.contains(x)))
            {
                phasedCounters.forEach(x -> ReadCounters.add(x));
                ++ReadCount;
                return true;
            }

            return false;
        }
    }
}
