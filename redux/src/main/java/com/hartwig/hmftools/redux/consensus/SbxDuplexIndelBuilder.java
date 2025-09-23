package com.hartwig.hmftools.redux.consensus;

import static java.lang.Math.abs;
import static java.lang.Math.min;
import static java.lang.String.format;

import static htsjdk.samtools.CigarOperator.I;
import static htsjdk.samtools.CigarOperator.S;

import java.util.Collections;
import java.util.List;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.bam.CigarUtils;
import com.hartwig.hmftools.common.bam.SupplementaryReadData;

import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMRecord;

public class SbxDuplexIndelBuilder
{
    private final SAMRecord mRecord;
    private final int mReadLength;
    private final List<CigarElement> mReadCigarElements;
    private final List<SbxDuplexIndel> mDuplexIndels;

    // use of supplementary data aligned soft-clip bases
    private SupplementaryReadData mSuppData;
    private List<CigarElement> mSuppDataCigarElements;
    private List<CigarElement> mAdjustedCigarElements;
    private int mSuppDataIndexStart;
    private int mSuppDataIndexEnd;
    private boolean mSuppDataInLeftSoftClip;

    public SbxDuplexIndelBuilder(final SAMRecord record, final List<Integer> duplexIndelIndices)
    {
        mRecord = record;
        mReadLength = mRecord.getReadBases().length;
        mDuplexIndels = Lists.newArrayListWithExpectedSize(duplexIndelIndices.size());

        mReadCigarElements = record.getCigar().getCigarElements();
        mSuppData = SupplementaryReadData.extractAlignment(record);

        mAdjustedCigarElements = mReadCigarElements;
        mSuppDataCigarElements = null;
        mSuppDataIndexStart = -1;
        mSuppDataIndexEnd = -1;

        if(mSuppData != null)
        {
            buildAdjustedSupplementaryCigar();
        }

        for(int i = 0; i < duplexIndelIndices.size();)
        {
            int duplexIndelIndexStart = duplexIndelIndices.get(i);
            int duplexIndelIndexEnd = duplexIndelIndexStart;

            ReadBaseInfo readBaseInfo = getReadBaseInfo(duplexIndelIndexStart);
            ReadBaseInfo nextReadBaseInfo = null;

            for(int j = i + 1; j < duplexIndelIndices.size(); ++j)
            {
                int nextDuplexIndelIndex = duplexIndelIndices.get(j);

                if(nextDuplexIndelIndex > duplexIndelIndexEnd + 1)
                    break;

                if(nextReadBaseInfo == null)
                    nextReadBaseInfo = new ReadBaseInfo(readBaseInfo);

                moveNext(nextReadBaseInfo);

                if(!nextReadBaseInfo.isReadBase())
                    break;

                duplexIndelIndexEnd = nextDuplexIndelIndex;
            }

            processDuplexIndel(readBaseInfo, duplexIndelIndexStart, duplexIndelIndexEnd);

            i += duplexIndelIndexEnd - duplexIndelIndexStart + 1;
        }
    }

    public List<SbxDuplexIndel> duplexIndels() { return mDuplexIndels; }

    private void processDuplexIndel(final ReadBaseInfo readBaseInfo, int duplexIndelIndexStart, int duplexIndelIndexEnd)
    {
        int duplexMismatchLength = duplexIndelIndexEnd - duplexIndelIndexStart + 1;
        byte[] repeatBases = new byte[duplexMismatchLength];

        for(int i = 0; i < duplexMismatchLength; ++i)
        {
            repeatBases[i] = mRecord.getReadBases()[duplexIndelIndexStart + i];
        }

        // convert to 1-base homopolymer if all bases match
        if(repeatBases.length > 1)
        {
            boolean isSingleBase = true;
            for(int i = 1; i < repeatBases.length; ++i)
            {
                if(repeatBases[i] != repeatBases[0])
                {
                    isSingleBase = false;
                    break;
                }
            }

            if(isSingleBase)
            {
                repeatBases = new byte[] { repeatBases[0] };
            }
        }

        int repeatLength = repeatBases.length;

        // search backwards and within inserted bases to find the start of the repeat
        int repeatStrIndex = repeatLength - 1;
        int insertRepeatCount = 0;
        int firstReadInsertIndex = -1;

        while(readBaseInfo.Index >= 0)
        {
            movePrevious(readBaseInfo);

            if(readBaseInfo.CigarOp != I)
                break;

            if(repeatBases[repeatStrIndex] != readBaseInfo.Base)
                break;

            --repeatStrIndex;

            if(repeatStrIndex < 0)
            {
                ++insertRepeatCount;
                repeatStrIndex = repeatLength - 1;
                firstReadInsertIndex = readBaseInfo.Index;
            }
        }

        // the repeat must be either wholy contained within the inserted bases or continue on past the insert without a gap
        if(firstReadInsertIndex < 0 || insertRepeatCount == 0)
            return;

        int insertRepeatLength = insertRepeatCount * repeatLength;

        int totalRepeatBaseLength = insertRepeatLength + duplexMismatchLength;

        int trimLength = min(insertRepeatLength, duplexMismatchLength);

        int lowBaseQualCount = duplexMismatchLength;

        int trimmedRepeatBaseLength = totalRepeatBaseLength - trimLength;

        // ensure there are enough low-qual bases for symmetry around the repeat
        if(lowBaseQualCount < trimmedRepeatBaseLength && (duplexMismatchLength % 2) == 1)
            ++lowBaseQualCount;

        int trimCount = 0;
        int deletedIndelIndexStart =-1;
        int deletedIndelIndexEnd =-1;

        moveTo(readBaseInfo, firstReadInsertIndex);

        List<Integer> lowQualIndices = Lists.newArrayListWithExpectedSize(lowBaseQualCount);
        int remainingLowQualCount = lowBaseQualCount;
        int addedLowQualAtStartCount = 0;

        for(int i = firstReadInsertIndex; i <= duplexIndelIndexEnd; ++i)
        {
            if(readBaseInfo.CigarOp == I && trimCount < trimLength) // only trim inserted bases
            {
                trimCount++;

                if(deletedIndelIndexStart < 0)
                {
                    deletedIndelIndexStart = i;
                    deletedIndelIndexEnd = i;
                }
                else
                {
                    deletedIndelIndexEnd = i;
                }
            }
            else if(remainingLowQualCount > 0)
            {
                lowQualIndices.add(i);
                --remainingLowQualCount;
                ++addedLowQualAtStartCount;

                if(remainingLowQualCount > 0)
                {
                    // apply symmetrically
                    int otherLowQualIndex = duplexIndelIndexEnd - (addedLowQualAtStartCount - 1);
                    lowQualIndices.add(otherLowQualIndex);
                    --remainingLowQualCount;
                }
            }
        }

        SbxDuplexIndel duplexIndel = new SbxDuplexIndel(
                duplexIndelIndexStart, duplexIndelIndexEnd, new String(repeatBases), totalRepeatBaseLength,
                firstReadInsertIndex, lowQualIndices, deletedIndelIndexStart, deletedIndelIndexEnd);

        mDuplexIndels.add(duplexIndel);
    }

    private void buildAdjustedSupplementaryCigar()
    {
        mAdjustedCigarElements = Lists.newArrayList();

        int leftSoftClipLength = mReadCigarElements.get(0).getOperator() == S ? mReadCigarElements.get(0).getLength() : 0;
        int lastElement = mReadCigarElements.size() - 1;
        int rightSoftClipLength = mReadCigarElements.get(lastElement).getOperator() == S ? mReadCigarElements.get(lastElement).getLength() : 0;

        if(leftSoftClipLength > rightSoftClipLength)
        {
            mSuppDataInLeftSoftClip = true;
            mSuppDataIndexStart = 0;
            mSuppDataIndexEnd = mSuppDataIndexStart + leftSoftClipLength - 1;
        }
        else
        {
            mSuppDataInLeftSoftClip = false;
            mSuppDataIndexEnd = mReadLength - 1;
            mSuppDataIndexStart = mSuppDataIndexEnd - rightSoftClipLength + 1;
        }

        List<CigarElement> suppElements = CigarUtils.cigarElementsFromStr(mSuppData.Cigar);

        if(mSuppData.isForwardOrient() == mRecord.getReadNegativeStrandFlag())
            Collections.reverse(suppElements);

        if(leftSoftClipLength > rightSoftClipLength)
        {
            mSuppDataCigarElements = checkSupplementaryCigar(suppElements, leftSoftClipLength, true);
            mAdjustedCigarElements.addAll(mSuppDataCigarElements);
            mAdjustedCigarElements.addAll(mReadCigarElements.subList(1, lastElement + 1)); // remove soft-clip from original
        }
        else
        {
            mAdjustedCigarElements.addAll(mReadCigarElements.subList(0, lastElement));

            mSuppDataCigarElements = checkSupplementaryCigar(suppElements, rightSoftClipLength, false);
            mAdjustedCigarElements.addAll(mSuppDataCigarElements);
        }
    }

    private boolean useLeftSoftClipSuppData() { return mSuppData != null && mSuppDataInLeftSoftClip; }
    private boolean useRightSoftClipSuppData() { return mSuppData != null && !mSuppDataInLeftSoftClip; }

    protected static List<CigarElement> checkSupplementaryCigar(
            final List<CigarElement> suppCigarElements, int softClipLength, boolean isLeftClipped)
    {
        // the supplementary data's cigar may not have the same mirrored soft-clip length, in which case it needs to be made to match
        int cigarCount = suppCigarElements.size();
        CigarElement matchingSoftClip;
        List<CigarElement> newSuppElements;

        // left clipped refers to the read, so the supp data cigar's right side is being evaluated and trimmed
        boolean trimLeftSide = !isLeftClipped;

        if(trimLeftSide)
        {
            matchingSoftClip = suppCigarElements.get(0);
            newSuppElements = suppCigarElements.subList(1, cigarCount);
        }
        else
        {
            matchingSoftClip = suppCigarElements.get(cigarCount - 1);
            newSuppElements = suppCigarElements.subList(0, cigarCount - 1);
        }

        // remove delete elements since these won't correspond to soft-clipped bases
        newSuppElements = newSuppElements.stream().filter(x -> x.getOperator().consumesReadBases()).collect(Collectors.toList());

        int readBaseLength = newSuppElements.stream().filter(x -> x.getOperator().consumesReadBases()).mapToInt(x -> x.getLength()).sum();

        if(readBaseLength == softClipLength)
            return newSuppElements;

        int diff = readBaseLength - softClipLength;

        if(diff > 0) // needs to be trimmed
        {
            trimCigar(newSuppElements, diff, trimLeftSide);
        }
        else
        {
            // add back the required portion of soft-clipping
            CigarElement extraElement = new CigarElement(abs(diff), matchingSoftClip.getOperator());

            if(trimLeftSide)
                newSuppElements.add(0, extraElement);
            else
                newSuppElements.add(extraElement);
        }

        return newSuppElements;
    }

    protected static void trimCigar(final List<CigarElement> cigarElements, int length, boolean fromStart)
    {
        int remainingBases = length;

        int index = fromStart ? 0 : cigarElements.size() - 1;

        while(remainingBases > 0)
        {
            CigarElement element = cigarElements.get(index);

            if(element.getLength() > remainingBases)
            {
                cigarElements.set(index, new CigarElement(element.getLength() - remainingBases, element.getOperator()));
                return;
            }
            else
            {
                cigarElements.remove(index);
                remainingBases -= element.getLength();

                if(!fromStart)
                    --index;
            }
        }
    }

    private class ReadBaseInfo
    {
        public int Index;
        public byte Base;
        public int CigarIndex;
        public CigarOperator CigarOp;
        public int IndexInCigar;
        public boolean Valid;

        public ReadBaseInfo(int index, byte base, int cigarIndex, final CigarOperator cigarOperator, int indexInCigar)
        {
            Index = index;
            Base = base;
            CigarOp = cigarOperator;
            CigarIndex = cigarIndex;
            IndexInCigar = indexInCigar;
            Valid = true;
        }

        public ReadBaseInfo(final ReadBaseInfo other)
        {
            Index = other.Index;
            Base = other.Base;
            CigarIndex = other.CigarIndex;
            CigarOp = other.CigarOp;
            IndexInCigar = other.IndexInCigar;
            Valid = other.Valid;
        }

        public boolean isReadBase() { return CigarOp.consumesReadBases(); }
        public boolean isRefBase() { return CigarOp.consumesReferenceBases(); }

        public String toString()
        {
            return format("%d: %c cigar(%s %d / %d)", Index, (char)Base, CigarOp, CigarIndex, IndexInCigar);
        }
    }

    private ReadBaseInfo getReadBaseInfo(int targetReadIndex)
    {
        int currentCigarIndex = 0;

        int readIndex = 0;
        for(; currentCigarIndex < mAdjustedCigarElements.size(); ++currentCigarIndex)
        {
            CigarElement element = mAdjustedCigarElements.get(currentCigarIndex);

            if(!element.getOperator().consumesReadBases())
                continue;

            if(readIndex + element.getLength() > targetReadIndex)
            {
                int indexInCigar = targetReadIndex - readIndex;

                return new ReadBaseInfo(
                        targetReadIndex, mRecord.getReadBases()[targetReadIndex], currentCigarIndex, element.getOperator(), indexInCigar);
            }

            readIndex += element.getLength();
        }

        return null;
    }

    private void moveTo(final ReadBaseInfo readBaseInfo, int targetReadIndex)
    {
        if(targetReadIndex == readBaseInfo.Index)
            return;

        if(targetReadIndex > readBaseInfo.Index)
        {
            while(targetReadIndex > readBaseInfo.Index)
            {
                moveNext(readBaseInfo);
            }
        }
        else
        {
            while(targetReadIndex > readBaseInfo.Index)
            {
                movePrevious(readBaseInfo);
            }
        }
    }

    private void movePrevious(final ReadBaseInfo readBaseInfo)
    {
        --readBaseInfo.Index;

        if(readBaseInfo.Index < 0)
        {
            readBaseInfo.Valid = false;
            return;
        }

        readBaseInfo.Base = mRecord.getReadBases()[readBaseInfo.Index];

        CigarElement element = mAdjustedCigarElements.get(readBaseInfo.CigarIndex);
        --readBaseInfo.IndexInCigar;

        if(readBaseInfo.IndexInCigar < 0)
        {
            --readBaseInfo.CigarIndex;

            if(readBaseInfo.CigarIndex < 0)
            {
                readBaseInfo.Valid = false;
                return;
            }

            element = mAdjustedCigarElements.get(readBaseInfo.CigarIndex);
            readBaseInfo.CigarOp = element.getOperator();
            readBaseInfo.IndexInCigar = element.getLength() - 1;
        }
    }

    private void moveNext(final ReadBaseInfo readBaseInfo)
    {
        ++readBaseInfo.Index;

        if(readBaseInfo.Index >= mReadLength)
        {
            readBaseInfo.Valid = false;
            return;
        }

        readBaseInfo.Base = mRecord.getReadBases()[readBaseInfo.Index];

        CigarElement element = mAdjustedCigarElements.get(readBaseInfo.CigarIndex);
        ++readBaseInfo.IndexInCigar;

        if(readBaseInfo.IndexInCigar >= element.getLength())
        {
            if(readBaseInfo.CigarIndex >= mAdjustedCigarElements.size())
            {
                readBaseInfo.Valid = false;
                return;
            }

            ++readBaseInfo.CigarIndex;
            element = mAdjustedCigarElements.get(readBaseInfo.CigarIndex);
            readBaseInfo.CigarOp = element.getOperator();
            readBaseInfo.IndexInCigar = 0;
        }
    }
}
