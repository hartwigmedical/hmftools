package com.hartwig.hmftools.sage.evidence;

import static java.lang.Math.abs;
import static java.lang.Math.max;
import static java.lang.Math.min;
import static java.lang.String.format;

import static com.hartwig.hmftools.common.region.BaseRegion.positionsOverlap;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_END;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_START;
import static com.hartwig.hmftools.sage.SageCommon.SG_LOGGER;
import static com.hartwig.hmftools.sage.evidence.FragmentSyncOutcome.ORIG_READ_COORDS;
import static com.hartwig.hmftools.sage.evidence.FragmentSyncType.BASE_MISMATCH;
import static com.hartwig.hmftools.sage.evidence.FragmentSyncType.CIGAR_MISMATCH;
import static com.hartwig.hmftools.sage.evidence.FragmentSyncType.COMBINED;
import static com.hartwig.hmftools.sage.evidence.FragmentSyncType.EXCEPTION;
import static com.hartwig.hmftools.sage.evidence.FragmentSyncType.NO_OVERLAP;
import static com.hartwig.hmftools.sage.evidence.FragmentSyncType.NO_OVERLAP_CIGAR_DIFF;
import static com.hartwig.hmftools.sage.evidence.FragmentSyncUtils.buildSyncedRead;
import static com.hartwig.hmftools.sage.evidence.FragmentSyncUtils.getCombinedBaseAndQual;
import static com.hartwig.hmftools.sage.evidence.FragmentSyncUtils.ignoreCigarOperatorMismatch;
import static com.hartwig.hmftools.sage.evidence.FragmentSyncUtils.isDeleteOrSplit;
import static com.hartwig.hmftools.sage.evidence.FragmentSyncUtils.overlappingCigarDiffs;

import static htsjdk.samtools.CigarOperator.D;
import static htsjdk.samtools.CigarOperator.I;
import static htsjdk.samtools.CigarOperator.M;
import static htsjdk.samtools.CigarOperator.N;
import static htsjdk.samtools.CigarOperator.S;

import java.util.Map;

import com.google.common.collect.Maps;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordSetBuilder;

public class FragmentSync
{
    private final FragmentSyncReadHandler mReadHandler;
    private final Map<String,SAMRecord> mCachedReads;

    private final int[] mSyncCounts;

    public FragmentSync(final FragmentSyncReadHandler readHandler)
    {
        mReadHandler = readHandler;
        mCachedReads = Maps.newHashMap();
        mSyncCounts = new int[FragmentSyncType.values().length];
    }

    public void clear() { mCachedReads.clear(); }
    public final int[] getSynCounts() { return mSyncCounts; }

    public boolean handleOverlappingReads(final SAMRecord record)
    {
        final SAMRecord otherRecord = mCachedReads.get(record.getReadName());

        if(otherRecord != null)
        {
            try
            {
                FragmentSyncOutcome syncOutcome = formFragmentRead(otherRecord, record);
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
                    mReadHandler.processReadRecord(fragmentRecord, false);
                }
                else if(syncOutcome.SyncType.processSeparately())
                {
                    // process both reads if a consensus failed
                    mReadHandler.processReadRecord(otherRecord, false);
                    mReadHandler.processReadRecord(record, false);
                }
                else if(syncOutcome.SyncType == CIGAR_MISMATCH)
                {
                    int firstIndelLen = otherRecord.getCigar().getCigarElements().stream()
                            .filter(x -> x.getOperator().isIndel()).mapToInt(x -> x.getLength()).sum();

                    int secondIndelLen = record.getCigar().getCigarElements().stream()
                            .filter(x -> x.getOperator().isIndel()).mapToInt(x -> x.getLength()).sum();

                    if(secondIndelLen < firstIndelLen)
                        mReadHandler.processReadRecord(record, false);
                    else
                        mReadHandler.processReadRecord(otherRecord, false);
                }
                else
                {
                    // only the first record
                    mReadHandler.processReadRecord(otherRecord, false);
                }
            }
            catch(Exception e)
            {
                ++mSyncCounts[FragmentSyncType.EXCEPTION.ordinal()];

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

        // no cache for reads where the mate doesn't overlap
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

    public static FragmentSyncOutcome formFragmentRead(final SAMRecord first, final SAMRecord second)
    {
        // take the highest base qual base for any overlapping bases
        // widen the read to cover both reads
        // how to handle preserve INDELs?

        int firstPosStart = first.getAlignmentStart();
        int firstPosEnd = first.getAlignmentEnd();
        int firstLength = first.getReadLength();

        int secondPosStart = second.getAlignmentStart();
        int secondPosEnd = second.getAlignmentEnd();
        int secondLength = second.getReadLength();

        if(!positionsOverlap(firstPosStart, firstPosEnd, secondPosStart, secondPosEnd))
        {
            return new FragmentSyncOutcome(NO_OVERLAP);
        }

        Cigar firstCigar = first.getCigar();
        Cigar secondCigar = second.getCigar();

        final byte[] firstBaseQualities = first.getBaseQualities();
        final byte[] firstBases = first.getReadBases();

        final CigarBaseCounts firstBaseCounts = new CigarBaseCounts(firstCigar);
        final CigarBaseCounts secondBaseCounts = new CigarBaseCounts(secondCigar);

        final byte[] secondBaseQualities = second.getBaseQualities();
        final byte[] secondBases = second.getReadBases();

        // work out boundaries and lengths
        int firstEffectivePosStart = firstPosStart - firstBaseCounts.SoftClipStart;
        int secondEffectivePosStart = secondPosStart - secondBaseCounts.SoftClipStart;
        int firstEffectivePosEnd = firstPosEnd + firstBaseCounts.SoftClipEnd;
        int secondEffectivePosEnd = secondPosEnd + secondBaseCounts.SoftClipEnd;

        if(firstBaseCounts.AdjustedBases != secondBaseCounts.AdjustedBases
        && overlappingCigarDiffs(firstCigar, firstPosStart, secondCigar, secondPosStart))
        {
            return new FragmentSyncOutcome(CIGAR_MISMATCH);
        }

        int combinedEffectiveStart = min(firstEffectivePosStart, secondEffectivePosStart);
        int combinedEffectiveEnd = max(firstEffectivePosEnd, secondEffectivePosEnd);

        // truncate any fragment with an insert size less than the expected read length
        int fragmentLength = abs(first.getInferredInsertSize());

        int truncatedFragmentStart = 0;
        int truncatedFragmentEnd = combinedEffectiveEnd;

        if(firstBaseCounts.AlignedBases > fragmentLength || secondBaseCounts.AdjustedBases > fragmentLength)
        {
            int excessBases = max(firstBaseCounts.AlignedBases, secondBaseCounts.AdjustedBases) - fragmentLength;
            truncatedFragmentStart = max(firstEffectivePosStart, secondEffectivePosStart);

            if(truncatedFragmentStart - combinedEffectiveStart < excessBases)
            {
                int endBaseTrim = excessBases - (truncatedFragmentStart - combinedEffectiveStart);
                truncatedFragmentEnd -= endBaseTrim;
            }
        }

        int combinedLength = combinedEffectiveEnd - combinedEffectiveStart + 1 + firstBaseCounts.AdjustedBases;

        final byte[] combinedBaseQualities = new byte[combinedLength];
        final byte[] combinedBases = new byte[combinedLength];
        int combinedPosStart = min(firstPosStart, secondPosStart);
        int combinedPosEnd = max(firstPosEnd, secondPosEnd);

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
            if(currentPos > combinedEffectiveStart && !isDeleteOrSplit(combinedCigarOperator))
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
                if(!isDeleteOrSplit(firstElement.getOperator()))
                    ++firstReadIndex;

                if(firstReadIndex >= firstLength)
                {
                    firstElement = null;
                }
                else if(firstReadIndex > firstCigarElementReadIndex + firstElement.getLength() - 1
                || (isDeleteOrSplit(firstElement.getOperator()) && combinedCigarElementLength == firstElement.getLength()))
                {
                    // move to next
                    if(!isDeleteOrSplit(firstElement.getOperator()))
                        firstCigarElementReadIndex += firstElement.getLength();

                    firstElement = firstCigarIndex < firstCigar.getCigarElements().size() - 1 ?
                            firstCigar.getCigarElement(++firstCigarIndex) : null;
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
                if(!isDeleteOrSplit(secondElement.getOperator()))
                    ++secondReadIndex;

                if(secondReadIndex >= secondLength)
                {
                    secondElement = null;
                }
                else if(secondReadIndex > secondCigarElementReadIndex + secondElement.getLength() - 1
                || (isDeleteOrSplit(secondElement.getOperator()) && combinedCigarElementLength == secondElement.getLength()))
                {
                    if(!isDeleteOrSplit(secondElement.getOperator()))
                        secondCigarElementReadIndex += secondElement.getLength();

                    secondElement = secondCigarIndex < secondCigar.getCigarElements().size() - 1 ?
                            secondCigar.getCigarElement(++secondCigarIndex) : null;
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
                        return new FragmentSyncOutcome(CIGAR_MISMATCH);
                    }
                }
                else if(firstElement != null && firstElement.getOperator().isIndel())
                {
                    return new FragmentSyncOutcome(NO_OVERLAP_CIGAR_DIFF);
                }
                else if(secondElement != null && secondElement.getOperator().isIndel())
                {
                    return new FragmentSyncOutcome(NO_OVERLAP_CIGAR_DIFF);
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
                    return new FragmentSyncOutcome(EXCEPTION);
                }
            }

            ++combinedCigarElementLength;

            if(combinedCigarOperator == I)
            {
                --currentPos;
            }
            else if(isDeleteOrSplit(combinedCigarOperator))
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
                        return new FragmentSyncOutcome(BASE_MISMATCH);
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

        SAMRecord combinedRecord = buildSyncedRead(
                first, combinedPosStart, combinedBases, combinedBaseQualities, combinedCigar,
                format("%d;%d;%d;%d", firstPosStart, firstPosEnd, secondPosStart, secondPosEnd));

        return new FragmentSyncOutcome(combinedRecord, COMBINED);
    }
}
