package com.hartwig.hmftools.sage.sync;

import static java.lang.Math.abs;
import static java.lang.Math.max;
import static java.lang.Math.min;

import static com.hartwig.hmftools.common.qual.BaseQualAdjustment.BASE_QUAL_MINIMUM;
import static com.hartwig.hmftools.common.region.BaseRegion.positionWithin;
import static com.hartwig.hmftools.common.region.BaseRegion.positionsOverlap;
import static com.hartwig.hmftools.sage.sync.FragmentSyncType.BASE_MISMATCH;
import static com.hartwig.hmftools.sage.sync.FragmentSyncType.CIGAR_MISMATCH;
import static com.hartwig.hmftools.sage.sync.FragmentSyncType.COMBINED;
import static com.hartwig.hmftools.sage.sync.FragmentSyncType.INVERSION;
import static com.hartwig.hmftools.sage.sync.FragmentSyncType.NO_OVERLAP;

import static htsjdk.samtools.CigarOperator.I;
import static htsjdk.samtools.CigarOperator.M;
import static htsjdk.samtools.CigarOperator.S;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.qual.BaseQualAdjustment;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordSetBuilder;

public class CombinedSyncData
{
    // workings and components of the fragment creation
    private int mFragmentStart;
    private int mFragmentEnd;
    private int mCombinedEffectiveStart;
    private int mFirstEffectivePosStart;
    private int mSecondEffectivePosStart;

    private byte[] mCombinedBases;
    private byte[] mCombinedBaseQualities;
    private final List<CigarElement> mCombinedCigar;

    public CombinedSyncData()
    {
        mCombinedCigar = Lists.newArrayList();
        mCombinedBases = null;
        mCombinedBaseQualities = null;
        mFragmentStart = 0;
        mFragmentEnd = 0;
        mFirstEffectivePosStart = 0;
        mSecondEffectivePosStart = 0;
        mCombinedEffectiveStart = 0;
    }

    public static FragmentSyncOutcome formFragmentRead(final SAMRecord first, final SAMRecord second, final SAMFileHeader samHeader)
    {
        if(!positionsOverlap(first.getAlignmentStart(), first.getAlignmentEnd(), second.getAlignmentStart(), second.getAlignmentEnd()))
        {
            return new FragmentSyncOutcome(NO_OVERLAP);
        }

        if(first.getReadNegativeStrandFlag() == second.getReadNegativeStrandFlag())
        {
            return new FragmentSyncOutcome(INVERSION);
        }

        CombinedSyncData combinedSyncData = new CombinedSyncData();

        if(!combinedSyncData.buildCombinedCigar(first, second))
        {
            return new FragmentSyncOutcome(CIGAR_MISMATCH);
        }

        if(!combinedSyncData.setBaseAndQualities(first, second))
        {
            return new FragmentSyncOutcome(BASE_MISMATCH);
        }

        return new FragmentSyncOutcome(combinedSyncData.buildSyncedRead(first, second, samHeader), COMBINED);
    }

    private SAMRecord buildSyncedRead(final SAMRecord first, final SAMRecord second, final SAMFileHeader samHeader)
    {
        SAMRecordSetBuilder recordBuilder = new SAMRecordSetBuilder();
        recordBuilder.setUnmappedHasBasesAndQualities(false);

        if(samHeader != null)
            recordBuilder.setHeader(samHeader);
        else
            recordBuilder.setHeader(first.getHeader());

        SAMRecord combinedRecord = recordBuilder.addFrag(
                first.getReadName(),
                first.getReferenceIndex(),
                mFragmentStart,
                first.getReadNegativeStrandFlag(),
                false,
                new Cigar(mCombinedCigar).toString(), "", 1, false);

        combinedRecord.setReadBases(mCombinedBases);
        combinedRecord.setAlignmentStart(mFragmentStart);

        combinedRecord.setBaseQualities(mCombinedBaseQualities);
        combinedRecord.setMateAlignmentStart(mFragmentStart);
        combinedRecord.setMateReferenceIndex(first.getMateReferenceIndex());

        combinedRecord.setFlags(first.getFlags());

        combinedRecord.setMappingQuality(max(first.getMappingQuality(), second.getMappingQuality()));

        // no need to compute since both records have the same value and it remains unchanged
        combinedRecord.setInferredInsertSize(abs(first.getInferredInsertSize()));

        for(SAMRecord.SAMTagAndValue tagAndValue : first.getAttributes())
        {
            combinedRecord.setAttribute(tagAndValue.tag, tagAndValue.value);
        }

        return combinedRecord;
    }

    private boolean buildCombinedCigar(final SAMRecord first, final SAMRecord second)
    {
        int firstPosStart = first.getAlignmentStart();
        int firstPosEnd = first.getAlignmentEnd();

        int secondPosStart = second.getAlignmentStart();
        int secondPosEnd = second.getAlignmentEnd();

        CigarState firstCigar = new CigarState(first.getCigar().getCigarElements());
        CigarState secondCigar = new CigarState(second.getCigar().getCigarElements());

        // work out boundaries and lengths
        mFirstEffectivePosStart = firstPosStart - firstCigar.softClipStart();
        mSecondEffectivePosStart = secondPosStart - secondCigar.softClipStart();
        int firstEffectivePosEnd = firstPosEnd + firstCigar.softClipEnd();
        int secondEffectivePosEnd = secondPosEnd + secondCigar.softClipEnd();

        int combinedEffectiveEnd;

        // for short fragments, less than standard read length, cap fragment boundaries from the five-prime read
        if(!first.getReadNegativeStrandFlag())
        {
            // first is the 5' read at the start, second is 5' from the right/end
            mFragmentStart = firstPosStart;
            mFragmentEnd = secondPosEnd;

            mCombinedEffectiveStart = mFirstEffectivePosStart;
            combinedEffectiveEnd = secondEffectivePosEnd;
        }
        else
        {
            mFragmentStart = secondPosStart;
            mFragmentEnd = firstPosEnd;

            mCombinedEffectiveStart = mSecondEffectivePosStart;
            combinedEffectiveEnd = firstEffectivePosEnd;
        }

        // skip past any bases prior to the fragment effective start
        if(mFirstEffectivePosStart < mCombinedEffectiveStart)
            firstCigar.moveToPosition(mFirstEffectivePosStart, mCombinedEffectiveStart);
        else if(mSecondEffectivePosStart < mCombinedEffectiveStart)
            secondCigar.moveToPosition(mSecondEffectivePosStart, mCombinedEffectiveStart);

        // now walk through the 2 cigars, building a combined one and checking for any incompatibilities
        int combinedCigarElementLength = 0;
        CigarOperator combinedCigarOperator = null;

        boolean firstActive = false;
        boolean secondActive = false;

        for(int currentPos = mCombinedEffectiveStart; currentPos <= combinedEffectiveEnd; ++currentPos)
        {
            boolean firstCigarChange = false;
            boolean secondCigarChange = false;

            if(currentPos >= mFirstEffectivePosStart)
            {
                if(!firstActive)
                {
                    firstActive = true;
                    firstCigarChange = true;
                }
                else if(!firstCigar.exhausted())
                {
                    firstCigarChange = firstCigar.moveNext();
                }
            }

            if(currentPos >= mSecondEffectivePosStart)
            {
                if(!secondActive)
                {
                    secondActive = true;
                    secondCigarChange = true;
                }
                else if(!secondCigar.exhausted())
                {
                    secondCigarChange = secondCigar.moveNext();
                }
            }

            if(firstCigarChange || secondCigarChange)
            {
                // scenarios:
                // both change and they are incompatible (anything except M and S differences at start or end)
                // both change, they disagree and are different from combined cigar
                // only 1 changes, and is different from combined cigar
                CigarOperator firstOperator = firstActive ? firstCigar.currentOperator() : null;
                CigarOperator secondOperator = secondActive ? secondCigar.currentOperator() : null;

                CigarOperator newOperator = null;

                if(firstOperator != null && secondOperator != null)
                {
                    // check for compatibility of current element operators
                    if(firstOperator != secondOperator)
                    {
                        if(!ignoreCigarOperatorMismatch(firstOperator, secondOperator))
                            return false; // a mismatch

                        // always favour an aligned element if it is within the permitted 5' aligned region
                        if(positionWithin(currentPos, mFragmentStart, mFragmentEnd))
                            newOperator = M;
                        else
                            newOperator = S;
                    }
                    else
                    {
                        newOperator = firstOperator;
                    }
                }
                else if(firstOperator != null)
                {
                    newOperator = firstOperator;
                }
                else // secondCigarChange
                {
                    newOperator = secondOperator;
                }

                if(combinedCigarElementLength > 0 && combinedCigarOperator != newOperator)
                {
                    mCombinedCigar.add(new CigarElement(combinedCigarElementLength, combinedCigarOperator));
                    combinedCigarElementLength = 0;
                }

                combinedCigarOperator = newOperator;
            }

            ++combinedCigarElementLength;

            if(combinedCigarOperator == I)
            {
                --currentPos;
            }
        }

        // add the last cigar element
        mCombinedCigar.add(new CigarElement(combinedCigarElementLength, combinedCigarOperator));

        return true;
    }

    private boolean setBaseAndQualities(final SAMRecord first, final SAMRecord second)
    {
        int combinedBaseLength = mCombinedCigar.stream().mapToInt(x -> x.getOperator().consumesReadBases() ? x.getLength() : 0).sum();

        mCombinedBases = new byte[combinedBaseLength];
        mCombinedBaseQualities = new byte[combinedBaseLength];
        
        int firstReadIndex = findFragmentReadIndexStart(first, mFirstEffectivePosStart);
        int secondReadIndex = findFragmentReadIndexStart(second, mSecondEffectivePosStart);

        final byte[] firstBaseQualities = first.getBaseQualities();
        final byte[] firstBases = first.getReadBases();
        final int firstLength = first.getReadLength();

        final byte[] secondBaseQualities = second.getBaseQualities();
        final byte[] secondBases = second.getReadBases();
        final int secondLength = second.getReadLength();

        CigarState combinedCigar = new CigarState(mCombinedCigar);
        int currentPosition = mCombinedEffectiveStart;

        for(int combinedReadIndex = 0; combinedReadIndex < combinedBaseLength;)
        {
            boolean moveReadBases = combinedCigar.currentOperator().consumesReadBases();
            boolean moveRefPosition = moveEffectivePosition(combinedCigar.currentOperator());

            if(moveReadBases)
            {
                if(firstReadIndex >= 0 && firstReadIndex < firstLength && secondReadIndex >= 0 && secondReadIndex < secondLength)
                {
                    if(firstBases[firstReadIndex] == secondBases[secondReadIndex])
                    {
                        mCombinedBases[combinedReadIndex] = firstBases[firstReadIndex];
                        mCombinedBaseQualities[combinedReadIndex] =
                                (byte)max(firstBaseQualities[firstReadIndex], secondBaseQualities[secondReadIndex]);
                    }
                    else
                    {
                        byte[] baseAndQual = getCombinedBaseAndQual(
                                firstBases[firstReadIndex], firstBaseQualities[firstReadIndex],
                                secondBases[secondReadIndex], secondBaseQualities[secondReadIndex]);

                        mCombinedBases[combinedReadIndex] = baseAndQual[0];
                        mCombinedBaseQualities[combinedReadIndex] = BaseQualAdjustment.adjustBaseQual(baseAndQual[1]);
                    }
                }
                else if(firstReadIndex >= 0 && firstReadIndex < firstLength)
                {
                    mCombinedBases[combinedReadIndex] = firstBases[firstReadIndex];
                    mCombinedBaseQualities[combinedReadIndex] = firstBaseQualities[firstReadIndex];
                }
                else if(secondReadIndex >= 0 && secondReadIndex < secondLength)
                {
                    mCombinedBases[combinedReadIndex] = secondBases[secondReadIndex];
                    mCombinedBaseQualities[combinedReadIndex] = secondBaseQualities[secondReadIndex];
                }
            }

            if(moveRefPosition)
                ++currentPosition;

            combinedCigar.moveNext();

            if(moveReadBases)
                ++combinedReadIndex;

            if(firstReadIndex < 0)
            {
                if(currentPosition == mFirstEffectivePosStart)
                    firstReadIndex = 0;
            }
            else if(moveReadBases)
            {
                ++firstReadIndex;
            }

            if(secondReadIndex < 0)
            {
                if(currentPosition == mSecondEffectivePosStart)
                    secondReadIndex = 0;
            }
            else if(moveReadBases)
            {
                ++secondReadIndex;
            }
        }

        return true;
    }

    private static boolean moveEffectivePosition(final CigarOperator operator)
    {
        return operator.consumesReferenceBases() || operator == S;
    }

    private int findFragmentReadIndexStart(final SAMRecord record, final int effectiveStart)
    {
        if(effectiveStart > mCombinedEffectiveStart)
            return -1;

        if(effectiveStart == mCombinedEffectiveStart)
            return 0;

        int readPosition = effectiveStart;
        int readIndex = 0;

        for(CigarElement element : record.getCigar().getCigarElements())
        {
            int positionMove = 0;

            if(moveEffectivePosition(element.getOperator()))
            {
                positionMove = min(mCombinedEffectiveStart - readPosition, element.getLength());

                readPosition += positionMove;
            }

            if(element.getOperator().consumesReadBases())
            {
                if(positionMove == 0)
                    readIndex += element.getLength();
                else
                    readIndex += positionMove;
            }

            if(readPosition >= mCombinedEffectiveStart)
                break;
        }

        return readIndex;
    }

    protected static boolean ignoreCigarOperatorMismatch(final CigarOperator first, final CigarOperator second)
    {
        return (first == M || first == S) && (second == M || second == S);
    }

    protected static byte[] getCombinedBaseAndQual(byte firstBase, byte firstQual, byte secondBase, byte secondQual)
    {
        if(firstBase == secondBase)
        {
            byte qual = (byte)max(firstQual, secondQual);
            return new byte[] { firstBase, qual };
        }
        else if(firstQual > secondQual)
        {
            // use the difference in quals
            return new byte[] { firstBase, (byte)max(BASE_QUAL_MINIMUM, firstQual - secondQual) };
        }
        else
        {
            return new byte[] { secondBase, (byte)max(BASE_QUAL_MINIMUM, secondQual - firstQual) };
        }
    }
}
