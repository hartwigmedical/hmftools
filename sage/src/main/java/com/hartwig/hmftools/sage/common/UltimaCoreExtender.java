package com.hartwig.hmftools.sage.common;

import static java.lang.Math.max;

import static htsjdk.samtools.CigarOperator.S;

import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;

import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMRecord;

public class UltimaCoreExtender
{
    private static final byte MISSING_BASE = -1;

    public static class UltimaCoreInfo
    {
        public final int ReadCoreStart;
        public final int ReadCoreEnd;
        public final ReadCigarInfo CigarInfo;

        public UltimaCoreInfo(final int readCoreStart, final int readCoreEnd, final ReadCigarInfo cigarInfo)
        {
            ReadCoreStart = readCoreStart;
            ReadCoreEnd = readCoreEnd;
            CigarInfo = cigarInfo;
        }
    }

    private static class AlignedBase
    {
        public final int RefPos;
        public final int ReadIndex;
        public final byte RefBase;
        public final byte ReadBase;
        public final CigarOperator CigarOp;

        public AlignedBase(int refPos, int readIndex, byte refBase, byte readBase, final CigarOperator cigarOp)
        {
            RefPos = refPos;
            ReadIndex = readIndex;
            RefBase = refBase;
            ReadBase = readBase;
            CigarOp = cigarOp;
        }

        // TODO: Improve this.
        @Override
        public String toString()
        {
            return "(" + RefPos + ", " + (char) RefBase + ", " + ReadIndex + ", " + (char) ReadBase + "," + CigarOp.toString() + ")";
        }
    }

    // TODO: Can probably simplify.
    public static UltimaCoreInfo extendCore(final SAMRecord read, final RefSequence refSequence, final int readAlignmentStart, final List<CigarElement> cigarElements, final int readCoreStart, final int readCoreEnd, final ReadCigarInfo readCigarInfo, final int flankSize)
    {
        final byte[] readBases = read.getReadBases();
        final int refStartIndex = readAlignmentStart - refSequence.Start;
        final byte[] refBases = refSequence.Bases;

        // align read bases against ref
        int readIndex = 0;
        int refIndex = refStartIndex;
        int pos = readAlignmentStart;
        boolean nonClippingOpSeen = false;
        final List<AlignedBase> alignedBases = Lists.newArrayList();
        for(int i = 0; i < cigarElements.size(); i++)
        {
            CigarElement cigarElement = cigarElements.get(i);
            boolean isRef = cigarElement.getOperator().consumesReferenceBases();
            boolean isRead = cigarElement.getOperator().consumesReadBases();
            if(nonClippingOpSeen)
            {
                if(cigarElement.getOperator().isClipping())
                {
                    break;
                }
            }
            else if(!cigarElement.getOperator().isClipping())
            {
                nonClippingOpSeen = true;
            }
            else if(cigarElement.getOperator() == S)
            {
                // move ref indices back
                refIndex -= cigarElement.getLength();
                if(refIndex < 0)
                {
                    return null;
                }

                pos -= cigarElement.getLength();
                isRef = true;
            }

            if(!isRef && !isRead)
            {
                continue;
            }

            if(isRef && isRead)
            {
                for(int j = 0; j < cigarElement.getLength(); j++)
                {
                    AlignedBase alignedBase = new AlignedBase(pos, readIndex, refBases[refIndex], readBases[readIndex],  cigarElement.getOperator());
                    alignedBases.add(alignedBase);
                    pos++;
                    refIndex++;
                    readIndex++;
                }

                continue;
            }

            if(isRef)
            {
                for(int j = 0; j < cigarElement.getLength(); j++)
                {
                    AlignedBase alignedBase = new AlignedBase(pos, readIndex - 1, refBases[refIndex], MISSING_BASE, cigarElement.getOperator());
                    alignedBases.add(alignedBase);
                    pos++;
                    refIndex++;
                }

                continue;
            }

            for(int j = 0; j < cigarElement.getLength(); j++)
            {
                AlignedBase alignedBase = new AlignedBase(pos - 1, readIndex, MISSING_BASE, readBases[readIndex], cigarElement.getOperator());
                alignedBases.add(alignedBase);
                readIndex++;
            }
        }

        // create dictionaries of refPos and readIndex to aligned base index
        Map<Integer, Integer> lookupFromRefPos = Maps.newHashMap();
        Map<Integer, Integer> lookupFromReadIndex = Maps.newHashMap();
        for(int i = 0; i < alignedBases.size(); i++)
        {
            AlignedBase alignedBase = alignedBases.get(i);
            if(alignedBase.RefBase != MISSING_BASE)
            {
                lookupFromRefPos.put(alignedBase.RefPos, i);
            }

            if(alignedBase.ReadBase != MISSING_BASE)
            {
                lookupFromReadIndex.put(alignedBase.ReadIndex, i);
            }
        }

        // get start and end index for the core
        Integer coreStart = lookupFromRefPos.get(readCigarInfo.CorePositionStart);
        Integer coreEnd = lookupFromRefPos.get(readCigarInfo.CorePositionEnd);
        if(coreStart == null || coreEnd == null)
        {
            // can't check if expansion is needed
            return null;
        }

        // extend so we are not cutting off homopolymers
        boolean leftExpanded = false;
        boolean rightExpanded = false;
        while(true)
        {
            AlignedBase coreStartBase = alignedBases.get(coreStart);
            AlignedBase coreEndBase = alignedBases.get(coreEnd);

            int firstRefBaseIndex = coreStartBase.RefBase == MISSING_BASE ? lookupFromRefPos.get(coreStartBase.RefPos + 1) : coreStart;
            int lastRefBaseIndex = lookupFromRefPos.get(coreEndBase.RefPos);

            int firstReadBaseIndex = coreStartBase.ReadBase == MISSING_BASE ? lookupFromReadIndex.get(coreStartBase.ReadIndex + 1) : coreStart;
            int lastReadBaseIndex = lookupFromReadIndex.get(coreEndBase.ReadIndex);

            Integer prevRefBaseIndex = lookupFromRefPos.get(alignedBases.get(firstRefBaseIndex).RefPos - 1);
            Integer nextRefBaseIndex = lookupFromRefPos.get(alignedBases.get(lastRefBaseIndex).RefPos + 1);

            Integer prevReadBaseIndex = lookupFromReadIndex.get(alignedBases.get(firstReadBaseIndex).ReadIndex - 1);
            Integer nextReadBaseIndex = lookupFromReadIndex.get(alignedBases.get(lastReadBaseIndex).ReadIndex + 1);

            if(prevRefBaseIndex == null || nextRefBaseIndex == null || prevReadBaseIndex == null || nextReadBaseIndex == null)
            {
                return null;
            }

            boolean expandLeft = alignedBases.get(firstRefBaseIndex).RefBase == alignedBases.get(prevRefBaseIndex).RefBase
                    || alignedBases.get(firstReadBaseIndex).ReadBase == alignedBases.get(prevReadBaseIndex).ReadBase;
            boolean expandRight = alignedBases.get(lastRefBaseIndex).RefBase == alignedBases.get(nextRefBaseIndex).RefBase
                    || alignedBases.get(lastReadBaseIndex).ReadBase == alignedBases.get(nextReadBaseIndex).ReadBase;
            if(!expandLeft && !expandRight)
            {
                break;
            }

            if(expandLeft)
            {
                leftExpanded = true;
                coreStart--;
                if(coreStart < 0)
                {
                    return null;
                }
            }

            if(expandRight)
            {
                rightExpanded = true;
                coreEnd++;
                if(coreEnd >= alignedBases.size())
                {
                    return null;
                }
            }
        }

        // go left until we have one ref and read base padding
        if(leftExpanded)
        {
            boolean refPaddingAdded = false;
            boolean readPaddingAdded = false;
            while(coreStart >= 1 && (!refPaddingAdded || !readPaddingAdded))
            {
                coreStart--;
                refPaddingAdded = refPaddingAdded || alignedBases.get(coreStart).RefBase != MISSING_BASE;
                readPaddingAdded = readPaddingAdded || alignedBases.get(coreStart).ReadBase != MISSING_BASE;
            }

            if(!refPaddingAdded || !readPaddingAdded)
            {
                return null;
            }
        }

        if(rightExpanded)
        {
            boolean refPaddingAdded = false;
            boolean readPaddingAdded = false;
            while(coreEnd < alignedBases.size() - 1 && (!refPaddingAdded || !readPaddingAdded))
            {
                coreEnd++;
                refPaddingAdded = refPaddingAdded || alignedBases.get(coreEnd).RefBase != MISSING_BASE;
                readPaddingAdded = readPaddingAdded || alignedBases.get(coreEnd).ReadBase != MISSING_BASE;
            }

            if(!refPaddingAdded || !readPaddingAdded)
            {
                return null;
            }
        }

        // find beginning of left flank
        int flankStart = coreStart - 1;
        int readBasesSeen = 0;
        while(flankStart >= 0)
        {
            if(alignedBases.get(flankStart).ReadBase != MISSING_BASE)
            {
                readBasesSeen++;
            }

            if(readBasesSeen == flankSize)
            {
                break;
            }

            flankStart--;
        }

        if(readBasesSeen < flankSize)
        {
            return null;
        }

        // push alignment of left flank out if the first flank base is in an indel
        while(flankStart >= 0 && (alignedBases.get(flankStart).RefBase == MISSING_BASE
                || alignedBases.get(flankStart).ReadBase == MISSING_BASE))
        {
            flankStart--;
        }

        if(flankStart < 0)
        {
            return null;
        }

        // find ending of right flank
        int flankEnd = coreEnd + 1;
        readBasesSeen = 0;
        while(flankEnd < alignedBases.size())
        {
            if(alignedBases.get(flankEnd).ReadBase != MISSING_BASE)
            {
                readBasesSeen++;
            }

            if(readBasesSeen == flankSize)
            {
                break;
            }

            flankEnd++;
        }

        if(readBasesSeen < flankSize)
        {
            return null;
        }

        // push alignment of right flank out if the last flank base is in an indel
        while(flankEnd < alignedBases.size() && (alignedBases.get(flankEnd).RefBase == MISSING_BASE
                || alignedBases.get(flankEnd).ReadBase == MISSING_BASE))
        {
            flankEnd++;
        }

        if(flankEnd >= alignedBases.size())
        {
            return null;
        }

        // get cigar of the flanks and core
        List<CigarElement> readCigarElements = Lists.newArrayList();
        CigarOperator currentOp = null;
        int currentLength = 0;
        for(int i = flankStart; i <= flankEnd; i++)
        {
            CigarOperator op = alignedBases.get(i).CigarOp;
            if(currentOp == null)
            {
                currentOp = op;
                currentLength = 1;
                continue;
            }

            if(op == currentOp)
            {
                currentLength++;
                continue;
            }

            readCigarElements.add(new CigarElement(currentLength, currentOp));
            currentOp = op;
            currentLength = 1;
        }

        if(currentLength >= 1)
        {
            readCigarElements.add(new CigarElement(currentLength, currentOp));
        }

        int newReadCoreStart = alignedBases.get(coreStart).ReadIndex + (alignedBases.get(coreStart).ReadBase == MISSING_BASE ? 1 : 0);
        int newReadCoreEnd = alignedBases.get(coreEnd).ReadIndex;
        int corePosStart = alignedBases.get(coreStart).RefPos + (alignedBases.get(coreStart).RefBase == MISSING_BASE ? 1 : 0);
        int corePosEnd = alignedBases.get(coreEnd).RefPos;
        return new UltimaCoreInfo(
                newReadCoreStart,
                newReadCoreEnd,
                new ReadCigarInfo(
                        readAlignmentStart,
                        readCigarElements,
                        max(alignedBases.get(flankStart).RefPos, readAlignmentStart),
                        alignedBases.get(flankEnd).RefPos,
                        corePosStart,
                        corePosEnd,
                        alignedBases.get(flankStart).ReadIndex,
                        alignedBases.get(flankEnd).ReadIndex));
    }
}