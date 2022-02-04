package com.hartwig.hmftools.sage.evidence;

import static com.hartwig.hmftools.common.utils.sv.BaseRegion.positionWithin;
import static com.hartwig.hmftools.sage.evidence.ReadMatchType.NO_SUPPORT;
import static com.hartwig.hmftools.sage.evidence.ReadMatchType.SUPPORT;

import java.io.File;
import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.utils.PerformanceCounter;
import com.hartwig.hmftools.common.utils.sv.ChrBaseRegion;
import com.hartwig.hmftools.sage.candidate.Candidate;
import com.hartwig.hmftools.sage.config.SageConfig;
import com.hartwig.hmftools.sage.phase.PhaseSetCounter;
import com.hartwig.hmftools.sage.quality.QualityCalculator;
import com.hartwig.hmftools.sage.quality.QualityRecalibrationMap;
import com.hartwig.hmftools.sage.common.RefSequence;
import com.hartwig.hmftools.sage.common.SamSlicer;
import com.hartwig.hmftools.sage.read.NumberEvents;

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

    // state per slice region
    private RefSequence mRefSequence;
    private QualityCalculator mQualityCalculator;
    private List<ReadContextCounter> mReadCounters;
    private int mLastCandidateIndex;

    private final VariantPhaser mVariantPhaser;

    private final List<PerformanceCounter> mPerfCounters;

    public static final int PC_PHASE_READS = 0;
    public static final int PC_FORM_LPS = 1;

    public ReadContextEvidence(
            final SageConfig config, final ReferenceSequenceFile refGenome,
            final Map<String,QualityRecalibrationMap> qualityRecalibrationMap, final PhaseSetCounter phaseSetCounter)
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
        mVariantPhaser = new VariantPhaser(phaseSetCounter);

        mPerfCounters = Lists.newArrayList();
        mPerfCounters.add(new PerformanceCounter("PhaseReads"));
        mPerfCounters.add(new PerformanceCounter("FormLPS"));
    }

    public List<PerformanceCounter> getPerfCounters() { return mPerfCounters; }

    public List<ReadContextCounter> collectEvidence(
            final List<Candidate> candidates, final String sample, final String bam, boolean checkPhasing)
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

        mVariantPhaser.initialise(sliceRegion, checkPhasing, mSageConfig.LogLpsData);

        mRefSequence = new RefSequence(sliceRegion, mRefGenome);

        QualityRecalibrationMap qrMap = mQualityRecalibrationMap.get(sample);
        mQualityCalculator = new QualityCalculator(mSageConfig.Quality, qrMap, mRefSequence.IndexedBases);

        final SamReader tumorReader = SamReaderFactory.makeDefault().validationStringency(mSageConfig.Stringency)
                .referenceSource(new ReferenceSource(mRefGenome)).open(new File(bam));

        final SamSlicer slicer = new SamSlicer(0, sliceRegion);

        mPerfCounters.get(PC_PHASE_READS).start();
        mPerfCounters.get(PC_PHASE_READS).resume();

        slicer.slice(tumorReader, this::processReadRecord);

        mPerfCounters.get(PC_PHASE_READS).stop();

        // assign local phase set IDs to all phased variants
        mPerfCounters.get(PC_FORM_LPS).start();
        mVariantPhaser.assignLocalPhaseSets();
        mPerfCounters.get(PC_FORM_LPS).stop();

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

        while(prevIndex >= 0 && !mReadCounters.isEmpty())
        {
            ReadContextCounter readCounter = mReadCounters.get(prevIndex);

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

        List<ReadContextCounter> posPhasedCounters = mVariantPhaser.enabled() ? Lists.newArrayList() : null;
        List<ReadContextCounter> negPhasedCounters = mVariantPhaser.enabled() ? Lists.newArrayList() : null;

        for(ReadContextCounter readCounter : readCounters)
        {
            ReadMatchType matchType = readCounter.processRead(record, mSageConfig, mQualityCalculator, numberOfEvents);

            if(mVariantPhaser.enabled())
            {
                if(matchType == SUPPORT)
                    posPhasedCounters.add(readCounter);
                else if(matchType == NO_SUPPORT)
                    negPhasedCounters.add(readCounter);
            }
        }

        mPerfCounters.get(PC_PHASE_READS).resume();
        mVariantPhaser.registeredPhasedVariants(posPhasedCounters, negPhasedCounters);
        mPerfCounters.get(PC_PHASE_READS).pause();
    }

}
