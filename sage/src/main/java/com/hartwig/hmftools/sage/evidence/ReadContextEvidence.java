package com.hartwig.hmftools.sage.evidence;

import static java.lang.Math.max;
import static java.lang.Math.min;
import static java.lang.String.format;

import static com.hartwig.hmftools.common.utils.sv.BaseRegion.positionWithin;
import static com.hartwig.hmftools.common.utils.sv.BaseRegion.positionsOverlap;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_END;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_START;
import static com.hartwig.hmftools.sage.SageCommon.SG_LOGGER;
import static com.hartwig.hmftools.sage.evidence.ReadMatchType.NO_SUPPORT;
import static com.hartwig.hmftools.sage.evidence.ReadMatchType.SUPPORT;
import static com.hartwig.hmftools.sage.evidence.SyncFragmentType.*;

import static htsjdk.samtools.CigarOperator.D;
import static htsjdk.samtools.CigarOperator.I;
import static htsjdk.samtools.CigarOperator.M;
import static htsjdk.samtools.CigarOperator.S;

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
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordSetBuilder;

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

    private final int[] mSyncCounts;

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
        mSyncCounts = new int[SyncFragmentType.values().length];
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

    public final int[] getSynCounts() { return mSyncCounts; }

    private boolean handleOverlappingReads(final SAMRecord record)
    {
        final SAMRecord otherRecord = mCachedReads.get(record.getReadName());

        if(otherRecord != null)
        {
            try
            {
                SyncFragmentOutcome syncOutcome = formFragmentRead(otherRecord, record);
                mCachedReads.remove(record.getReadName());

                SAMRecord fragmentRecord = syncOutcome.CombinedRecord;
                ++mSyncCounts[syncOutcome.SyncType.ordinal()];

                /*
                SG_LOGGER.trace("fragment sync: first({} {}:{}-{} {}) second({} {}:{}-{} {}) outcome({}) {}",
                        otherRecord.getReadName(), otherRecord.getContig(), otherRecord.getAlignmentStart(), otherRecord.getAlignmentEnd(),
                        otherRecord.getCigarString(), record.getReadName(), record.getContig(), record.getAlignmentStart(),
                        record.getAlignmentEnd(), record.getCigarString(), syncOutcome.SyncType,
                        syncOutcome.CombinedRecord != null ? format("newCigar(%s)", syncOutcome.CombinedRecord.getCigarString()) : "");
                */

                if(fragmentRecord != null)
                {
                    processReadRecord(fragmentRecord, false);
                }
                else if(syncOutcome.SyncType.processSeparately())
                {
                    // process both reads if a consensus failed
                    processReadRecord(otherRecord, false);
                    processReadRecord(record, false);
                }
                else if(syncOutcome.SyncType == CIGAR_MISMATCH)
                {
                    int firstIndelLen = otherRecord.getCigar().getCigarElements().stream()
                            .filter(x -> x.getOperator().isIndel()).mapToInt(x -> x.getLength()).sum();

                    int secondIndelLen = record.getCigar().getCigarElements().stream()
                            .filter(x -> x.getOperator().isIndel()).mapToInt(x -> x.getLength()).sum();

                    if(secondIndelLen < firstIndelLen)
                        processReadRecord(record, false);
                    else
                        processReadRecord(otherRecord, false);
                }
                else
                {
                    // only the first record
                    processReadRecord(otherRecord, false);
                }
            }
            catch(Exception e)
            {
                ++mSyncCounts[SyncFragmentType.EXCEPTION.ordinal()];

                SG_LOGGER.error("failed to sync fragments: {}", e.toString());
                e.printStackTrace();

                SG_LOGGER.info("firstRead({} {}:{}-{} {})",
                        otherRecord.getReadName(), otherRecord.getContig(), otherRecord.getAlignmentStart(), otherRecord.getAlignmentEnd(),
                        otherRecord.getCigarString());

                SG_LOGGER.info("secondRead({} {}:{}-{} {})",
                        record.getReadName(), record.getContig(), record.getAlignmentStart(), record.getAlignmentEnd(),
                        record.getCigarString());
            }
        }

        if(!record.getContig().equals(record.getMateReferenceName()))
            return false;

        if(!positionsOverlap(
                record.getAlignmentStart(), record.getAlignmentEnd(),
                record.getMateAlignmentStart(), record.getMateAlignmentStart() + record.getReadLength()))
        {
            ++mSyncCounts[NO_OVERLAP.ordinal()];
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

        if(record.getCigar().getFirstCigarElement().getOperator() == S)
        {
            readStart -= record.getCigar().getFirstCigarElement().getLength();;
            readStart -= mMaxDeleteLength; // account for deleted bases being the cause of the soft-clipping
        }

        if(record.getCigar().getLastCigarElement().getOperator() == S)
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

    public static SyncFragmentOutcome formFragmentRead(final SAMRecord first, final SAMRecord second)
    {
        // take the highest base qual base for any overlapping bases
        // widen the read to cover both reads
        // how to handle preserve INDELs?

        int firstPosStart = first.getAlignmentStart();
        int firstPosEnd = first.getAlignmentEnd();
        int firstLength = first.getReadLength();

        int secondPosStart = second.getAlignmentStart();
        int secondPosEnd = second.getAlignmentEnd();
        int secondLength = first.getReadLength();

        if(!positionsOverlap(firstPosStart, firstPosEnd, secondPosStart, secondPosEnd))
        {
            return new SyncFragmentOutcome(NO_OVERLAP);
        }

        Cigar firstCigar = first.getCigar();
        Cigar secondCigar = second.getCigar();

        // must have matching non-alignment and SC elements, otherwise give up
        /*
        if(!compatibleCigars(firstCigar, secondCigar))
        {
            return new SyncFragmentOutcome(CIGAR_MISMATCH);
        }
        */

        final byte[] firstBaseQualities = first.getBaseQualities();
        final byte[] firstBases = first.getReadBases();

        int[] firstScLengths = new int[] {
                firstCigar.getFirstCigarElement().getOperator() == S ? firstCigar.getFirstCigarElement().getLength() : 0,
                firstCigar.getLastCigarElement().getOperator() == S ? firstCigar.getLastCigarElement().getLength() : 0
        };

        final byte[] secondBaseQualities = second.getBaseQualities();
        final byte[] secondBases = second.getReadBases();

        int[] secondScLengths = new int[] {
                secondCigar.getFirstCigarElement().getOperator() == S ? secondCigar.getFirstCigarElement().getLength() : 0,
                secondCigar.getLastCigarElement().getOperator() == S ? secondCigar.getLastCigarElement().getLength() : 0
        };

        // work out boundaries and lengths
        int firstEffectivePosStart = firstPosStart - firstScLengths[SE_START];
        int secondEffectivePosStart = secondPosStart - secondScLengths[SE_START];
        int firstEffectivePosEnd = firstPosEnd + firstScLengths[SE_END];
        int secondEffectivePosEnd = secondPosEnd + secondScLengths[SE_END];

        int combinedEffectiveStart = min(firstEffectivePosStart, secondEffectivePosStart);
        int combinedEffectiveEnd = max(firstEffectivePosEnd, secondEffectivePosEnd);

        int firstAdjustedBases = firstCigar.getCigarElements().stream()
                .filter(x -> x.getOperator().isIndel())
                .mapToInt(x -> x.getOperator() == D ? -x.getLength() : x.getLength())
                .sum();

        int secondAdjustedBases = secondCigar.getCigarElements().stream()
                .filter(x -> x.getOperator().isIndel())
                .mapToInt(x -> x.getOperator() == D ? -x.getLength() : x.getLength())
                .sum();

        if(firstAdjustedBases != secondAdjustedBases && overlappingCigarDiffs(firstCigar, firstPosStart, secondCigar, secondPosStart))
        {
            return new SyncFragmentOutcome(CIGAR_MISMATCH);
        }

        int combinedLength = combinedEffectiveEnd - combinedEffectiveStart + 1 + firstAdjustedBases;

        final byte[] combinedBaseQualities = new byte[combinedLength];
        final byte[] combinedBases = new byte[combinedLength];
        int combinedPosStart = min(firstPosStart, secondPosStart);

        int combinedReadIndex = 0;
        int firstReadIndex = -1;
        int secondReadIndex = -1;

        int firstCigarIndex = 0;
        CigarElement firstElement = null;
        int firstCigarElementReadIndex = 0;

        int secondCigarIndex = 0;
        CigarElement secondElement = null;
        int secondCigarElementReadIndex = 0;

        int combinedCigarElementLength = 0;
        CigarOperator combinedCigarOperator = M;
        Cigar combinedCigar = new Cigar();

        int baseMismatches = 0;

        for(int currentPos = combinedEffectiveStart; currentPos <= combinedEffectiveEnd; ++currentPos)
        {
            if(currentPos > combinedEffectiveStart && combinedCigarOperator != D)
                ++combinedReadIndex;

            boolean firstCigarChange = false;
            boolean secondCigarChange = false;

            if(currentPos == firstEffectivePosStart)
            {
                firstReadIndex = 0;
                firstElement = firstCigar.getCigarElement(firstCigarIndex);
                firstCigarChange = true;
            }
            else if(currentPos > firstEffectivePosStart && firstElement != null)
            {
                if(firstElement.getOperator() != D)
                    ++firstReadIndex;

                if(firstReadIndex >= firstLength)
                {
                    firstElement = null;
                }
                else if(firstReadIndex > firstCigarElementReadIndex + firstElement.getLength() - 1
                || (firstElement.getOperator() == D && combinedCigarElementLength == firstElement.getLength()))
                {
                    // move to next
                    if(firstElement.getOperator() != D)
                        firstCigarElementReadIndex += firstElement.getLength();

                    firstElement = firstCigar.getCigarElement(++firstCigarIndex);
                    firstCigarChange = true;
                }
            }

            if(currentPos == secondEffectivePosStart)
            {
                secondReadIndex = 0;
                secondElement = secondCigar.getCigarElement(secondCigarIndex);
                secondCigarChange = true;
            }
            else if(currentPos > secondEffectivePosStart && secondElement != null)
            {
                if(secondElement.getOperator() != D)
                    ++secondReadIndex;

                if(secondReadIndex >= secondLength)
                {
                    secondElement = null;
                }
                else if(secondReadIndex > secondCigarElementReadIndex + secondElement.getLength() - 1
                || (secondElement.getOperator() == D && combinedCigarElementLength == secondElement.getLength()))
                {
                    if(secondElement.getOperator() != D)
                        secondCigarElementReadIndex += secondElement.getLength();

                    secondElement = secondCigar.getCigarElement(++secondCigarIndex);
                    secondCigarChange = true;
                }
            }

            // more precise check of matching cigars
            if(firstCigarChange || secondCigarChange)
            {
                if(firstElement != null && secondElement != null)
                {
                    if(firstElement.getOperator() != secondElement.getOperator()
                    && !ignoreCigarOperatorMismatch(firstElement.getOperator(), secondElement.getOperator()))
                    {
                        return new SyncFragmentOutcome(CIGAR_MISMATCH);
                    }
                }
                else if(firstElement != null && firstElement.getOperator().isIndel())
                {
                    return new SyncFragmentOutcome(NO_OVERLAP_CIGAR_DIFF);
                }
                else if(secondElement != null && secondElement.getOperator().isIndel())
                {
                    return new SyncFragmentOutcome(NO_OVERLAP_CIGAR_DIFF);
                }
            }

            // handle cigar elements
            if((firstCigarChange || firstElement == null) && (secondCigarChange || secondElement == null))
            {
                if(combinedCigarElementLength > 0)
                {
                    combinedCigar.add(new CigarElement(combinedCigarElementLength, combinedCigarOperator));
                    combinedCigarElementLength = 0;
                }

                if(firstElement != null && secondElement != null)
                {
                    combinedCigarOperator = firstElement.getOperator() == M || secondElement.getOperator() == M ?
                            M : firstElement.getOperator();
                }
                else if(firstElement != null)
                {
                    combinedCigarOperator = firstElement.getOperator();
                }
                else if(secondElement != null)
                {
                    combinedCigarOperator = secondElement.getOperator();
                }
                else
                {
                    return new SyncFragmentOutcome(EXCEPTION);
                }
            }

            ++combinedCigarElementLength;

            if(combinedCigarOperator == I)
            {
                --currentPos;
            }
            else if(combinedCigarOperator == D)
            {
                if(combinedCigarElementLength <= firstElement.getLength())
                    continue;
            }

            // choose the base and qual
            if(firstReadIndex >= 0 && firstReadIndex < firstLength && secondReadIndex >= 0 && secondReadIndex < secondLength)
            {
                if(firstBases[firstReadIndex] == secondBases[secondReadIndex])
                {
                    combinedBases[combinedReadIndex] = firstBases[firstReadIndex];
                    combinedBaseQualities[combinedReadIndex] = (byte)max(firstBaseQualities[firstReadIndex], secondBaseQualities[secondReadIndex]);
                }
                else
                {
                    ++baseMismatches;

                    if(baseMismatches >= 10)
                    {
                        return new SyncFragmentOutcome(BASE_MISMATCH);
                    }

                    byte[] baseAndQual = getCombinedBaseAndQual(
                            firstBases[firstReadIndex], firstBaseQualities[firstReadIndex],
                            secondBases[secondReadIndex], secondBaseQualities[secondReadIndex]);

                    combinedBases[combinedReadIndex] = baseAndQual[0];
                    combinedBaseQualities[combinedReadIndex] = baseAndQual[1];
                }
            }
            else if(firstReadIndex >= 0 && firstReadIndex < firstLength)
            {
                combinedBases[combinedReadIndex] = firstBases[firstReadIndex];
                combinedBaseQualities[combinedReadIndex] = firstBaseQualities[firstReadIndex];
            }
            else
            {
                if(combinedReadIndex >= combinedBases.length || secondReadIndex >= secondBases.length)
                {
                    SG_LOGGER.error("here");
                    return null;
                }

                combinedBases[combinedReadIndex] = secondBases[secondReadIndex];
                combinedBaseQualities[combinedReadIndex] = secondBaseQualities[secondReadIndex];
            }
        }

        // add the last cigar element
        combinedCigar.add(new CigarElement(combinedCigarElementLength, combinedCigarOperator));

        SAMRecordSetBuilder recordBuilder = new SAMRecordSetBuilder();
        recordBuilder.setUnmappedHasBasesAndQualities(false);

        SAMRecord combinedRecord = recordBuilder.addFrag(
                first.getReadName(),
                first.getReferenceIndex(),
                combinedPosStart,
                first.getReadNegativeStrandFlag(),
                false,
                combinedCigar.toString(), "", 1, false);

        combinedRecord.setReadBases(combinedBases);
        combinedRecord.setAlignmentStart(combinedPosStart);
        combinedRecord.setReferenceIndex(first.getReferenceIndex());

        combinedRecord.setBaseQualities(combinedBaseQualities);
        combinedRecord.setReferenceName(first.getReferenceName());
        combinedRecord.setMateAlignmentStart(secondPosStart);
        combinedRecord.setMateReferenceName(second.getReferenceName());
        combinedRecord.setMateReferenceIndex(second.getReferenceIndex());

        // to be correct this should match the cigar element count
        combinedRecord.setFlags(first.getFlags());

        /*
        combinedRecord.setFirstOfPairFlag(true);
        combinedRecord.setReadPairedFlag(true);
        combinedRecord.setProperPairFlag(true);
        */

        combinedRecord.setMappingQuality(first.getMappingQuality());
        combinedRecord.setInferredInsertSize(combinedEffectiveEnd - combinedEffectiveStart + 1);

        for(SAMRecord.SAMTagAndValue tagAndValue : first.getAttributes())
        {
            combinedRecord.setAttribute(tagAndValue.tag, tagAndValue.value);
        }

        return new SyncFragmentOutcome(combinedRecord, COMBINED);
    }

    private static boolean ignoreCigarOperatorMismatch(final CigarOperator first, final CigarOperator second)
    {
        return (first == M || first == S) && (second == M || second == S);
    }

    public static boolean overlappingCigarDiffs(final Cigar firstCigar, int firstPosStart, final Cigar secondCigar, int secondPosStart)
    {
        int firstAdjustedElementPosEnd = 0;
        int readPos = firstPosStart;
        for(CigarElement element : firstCigar.getCigarElements())
        {
            switch(element.getOperator())
            {
                case M:
                    readPos += element.getLength();
                    break;
                case D:
                    readPos += element.getLength();
                    firstAdjustedElementPosEnd = readPos + element.getLength();
                    break;
                case I:
                    firstAdjustedElementPosEnd = readPos + 1;
                default:
                    break;
            }
        }

        int secondAdjustedElementPosStart = secondPosStart;
        for(CigarElement element : secondCigar.getCigarElements())
        {
            if(element.getOperator() == M)
                secondAdjustedElementPosStart += element.getLength();
            else if(element.getOperator().isIndel())
                break;
        }

        return firstAdjustedElementPosEnd >= secondAdjustedElementPosStart;
    }

    public static boolean compatibleCigars(final Cigar firstCigar, final Cigar secondCigar)
    {
        // SC at the start and end are optional, but otherwise all elements must match length and type
        int j = 0;
        int i = 0;

        CigarElement firstElement = firstCigar.getCigarElements().get(i);
        CigarElement secondElement = secondCigar.getCigarElements().get(j);

        if(firstElement.getOperator() == S)
            ++i;

        if(secondElement.getOperator() == S)
            ++j;

        while(true)
        {
            firstElement = i < firstCigar.getCigarElements().size() ? firstCigar.getCigarElements().get(i) : null;
            secondElement = j < secondCigar.getCigarElements().size() ? secondCigar.getCigarElements().get(j) : null;

            if(firstElement == null && secondElement == null)
                break;

            if(firstElement == null)
                return secondElement.getOperator() == S;
            else if(secondElement == null)
                return firstElement.getOperator() == S;

            // must match types and lengths if not an alignment
            if(firstElement.getOperator() != secondElement.getOperator())
                return false;

            if(firstElement.getOperator() == S)
                return true;

            if(firstElement.getOperator() != M && firstElement.getLength() != secondElement.getLength())
                return false;

            ++i;
            ++j;
        }

        return true;
    }

    public static byte[] getCombinedBaseAndQual(byte firstBase, byte firstQual, byte secondBase, byte secondQual)
    {
        if(firstBase == secondBase)
        {
            byte qual = (byte)max(firstQual, secondQual);
            return new byte[] { firstBase, qual };
        }
        else if(firstQual > secondQual)
        {
            // use the difference in quals
            return new byte[] { firstBase, (byte)((int)firstQual - (int)secondQual) };
        }
        else
        {
            return new byte[] { secondBase, (byte)((int)secondQual - (int)firstQual) };
        }
    }
}
