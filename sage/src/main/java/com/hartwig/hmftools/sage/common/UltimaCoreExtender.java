package com.hartwig.hmftools.sage.common;

import static java.lang.Math.max;
import static java.lang.String.format;

import static com.hartwig.hmftools.common.bam.CigarUtils.collapseCigarOps;

import static htsjdk.samtools.CigarOperator.M;
import static htsjdk.samtools.CigarOperator.S;

import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.utils.IntPair;

import org.jetbrains.annotations.Nullable;

import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;

public class UltimaCoreExtender
{
    @VisibleForTesting
    public static final byte MISSING_BASE = -1;
    @VisibleForTesting
    public static final int INVALID_INDEX = -1;

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

    @VisibleForTesting
    public static class AlignedBase
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

        @Override
        public String toString()
        {
            char refBase = RefBase == MISSING_BASE ? '_' : (char) RefBase;
            char readBase = ReadBase == MISSING_BASE ? '_' : (char) ReadBase;
            return format("ref(pos=%d base=%s) read(idx=%d base=%s) cigarOp(%s)", RefPos, refBase, ReadIndex, readBase, CigarOp.toString());
        }

        @Override
        public boolean equals(final Object o)
        {
            if(this == o)
                return true;

            if(!(o instanceof AlignedBase))
                return false;

            final AlignedBase that = (AlignedBase) o;
            return RefPos == that.RefPos && ReadIndex == that.ReadIndex && RefBase == that.RefBase && ReadBase == that.ReadBase
                    && CigarOp == that.CigarOp;
        }

        @Override
        public int hashCode()
        {
            int hash = RefPos;
            hash = 31 * hash + ReadIndex;
            hash = 31 * hash + RefBase;
            hash = 31 * hash + ReadBase;
            hash = 31 * hash + CigarOp.ordinal();
            return hash;
        }

        public boolean isExactMatch()
        {
            return RefBase != MISSING_BASE && RefBase == ReadBase && CigarOp == M;
        }

        public boolean isIndel()
        {
            return RefBase == MISSING_BASE || ReadBase == MISSING_BASE;
        }
    }

    @Nullable
    public static UltimaCoreInfo extendUltimaCore(final byte[] readBases, final RefSequence refSequence, final int readAlignmentStart,
            final List<CigarElement> cigarElements, final ReadCigarInfo readCigarInfo, final int flankSize, final boolean inAppendMode)
    {
        final List<AlignedBase> alignedBases = alignReadBases(readBases, refSequence, readAlignmentStart, cigarElements);
        if(alignedBases == null)
            return null;

        Map<Integer, Integer> lookupFromRefPos = Maps.newHashMap();
        Map<Integer, Integer> lookupFromReadIndex = Maps.newHashMap();
        populateAlignedBaseLookupMaps(alignedBases, lookupFromRefPos, lookupFromReadIndex);

        Integer coreStart = lookupFromRefPos.get(readCigarInfo.CorePositionStart);
        Integer coreEnd = lookupFromRefPos.get(readCigarInfo.CorePositionEnd);
        if(coreStart == null || coreEnd == null)
            return null;  // can't check if expansion is needed

        IntPair extendedReadCoreIndices = extendCore(alignedBases, coreStart, coreEnd, lookupFromRefPos, lookupFromReadIndex);
        if(extendedReadCoreIndices == null)
            return null;

        int newCoreStart = extendedReadCoreIndices.getLeft();
        int newCoreEnd = extendedReadCoreIndices.getRight();
        boolean addLeftPadding = newCoreStart < coreStart || !alignedBases.get(coreStart).isExactMatch();
        boolean addRightPadding = newCoreEnd > coreEnd || !alignedBases.get(coreEnd).isExactMatch();
        coreStart = newCoreStart;
        coreEnd = newCoreEnd;
        if(alignedBases.get(coreStart).CigarOp.isClipping() || alignedBases.get(coreEnd).CigarOp.isClipping())
            return null;

        if(addLeftPadding)
        {
            coreStart = addPadding(alignedBases, coreStart, true);
            if(coreStart == INVALID_INDEX)
                return null;

            if(alignedBases.get(coreStart).CigarOp.isClipping())
                return null;
        }

        if(addRightPadding)
        {
            coreEnd = addPadding(alignedBases, coreEnd, false);
            if(coreEnd == INVALID_INDEX)
                return null;

            if(alignedBases.get(coreEnd).CigarOp.isClipping())
                return null;
        }

        int flankStart = findFlankBoundary(alignedBases, coreStart, true, flankSize, inAppendMode);
        if(flankStart == INVALID_INDEX)
            return null;

        int flankEnd = findFlankBoundary(alignedBases, coreEnd, false, flankSize, inAppendMode);
        if(flankEnd == INVALID_INDEX)
            return null;

        // get cigar of the flanks and core
        List<CigarOperator> readCigarOps = IntStream.range(flankStart, flankEnd + 1)
                .mapToObj(i -> alignedBases.get(i).CigarOp)
                .collect(Collectors.toList());
        List<CigarElement> readCigarElements = collapseCigarOps(readCigarOps);

        int newReadCoreStart = alignedBases.get(coreStart).ReadIndex;
        int newReadCoreEnd = alignedBases.get(coreEnd).ReadIndex;
        int corePosStart = alignedBases.get(coreStart).RefPos;
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

    @VisibleForTesting
    public static void populateAlignedBaseLookupMaps(final List<AlignedBase> alignedBases, final Map<Integer, Integer> lookupFromRefPosOut,
            final Map<Integer, Integer> lookupFromReadIndexOut)
    {
        for(int i = 0; i < alignedBases.size(); i++)
        {
            AlignedBase alignedBase = alignedBases.get(i);
            if(alignedBase.RefBase != MISSING_BASE)
            {
                lookupFromRefPosOut.put(alignedBase.RefPos, i);
            }

            if(alignedBase.ReadBase != MISSING_BASE)
            {
                lookupFromReadIndexOut.put(alignedBase.ReadIndex, i);
            }
        }
    }

    @VisibleForTesting
    @Nullable
    public static List<AlignedBase> alignReadBases(final byte[] readBases, final RefSequence refSequence, final int readAlignmentStart,
            final List<CigarElement> cigarElements)
    {
        final byte[] refBases = refSequence.Bases;
        int readIndex = 0;
        int refIndex = readAlignmentStart - refSequence.Start;
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
                    break;
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
                    return null;

                pos -= cigarElement.getLength();
                isRef = true;
            }

            if(!isRef && !isRead)
                continue;

            if(isRef && isRead)
            {
                for(int j = 0; j < cigarElement.getLength(); j++)
                {
                    AlignedBase alignedBase = new AlignedBase(
                            pos, readIndex, refBases[refIndex], readBases[readIndex], cigarElement.getOperator());

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
                    AlignedBase alignedBase = new AlignedBase(
                            pos, readIndex - 1, refBases[refIndex], MISSING_BASE, cigarElement.getOperator());

                    alignedBases.add(alignedBase);
                    pos++;
                    refIndex++;
                }

                continue;
            }

            for(int j = 0; j < cigarElement.getLength(); j++)
            {
                AlignedBase alignedBase = new AlignedBase(
                        pos - 1, readIndex, MISSING_BASE, readBases[readIndex], cigarElement.getOperator());

                alignedBases.add(alignedBase);
                readIndex++;
            }
        }

        return alignedBases;
    }

    @VisibleForTesting
    @Nullable
    public static IntPair extendCore(final List<AlignedBase> alignedBases, final int initCoreStart, final int initCoreEnd,
            final Map<Integer, Integer> lookupFromRefPos, final Map<Integer, Integer> lookupFromReadIndex)
    {
        if(initCoreStart < 0 || initCoreEnd >= alignedBases.size())
            return null;

        int coreStart = initCoreStart;
        int coreEnd = initCoreEnd;
        while(true)
        {
            AlignedBase coreStartBase = alignedBases.get(coreStart);
            AlignedBase coreEndBase = alignedBases.get(coreEnd);

            if(coreStartBase.isIndel())
            {
                coreStart--;
                if(coreStart < 0)
                    return null;

                continue;
            }

            if(coreEndBase.isIndel())
            {
                coreEnd++;
                if(coreEnd >= alignedBases.size())
                    return null;

                continue;
            }

            Integer prevRefBaseIndex = lookupFromRefPos.get(coreStartBase.RefPos - 1);
            Integer nextRefBaseIndex = lookupFromRefPos.get(coreEndBase.RefPos + 1);

            Integer prevReadBaseIndex = lookupFromReadIndex.get(coreStartBase.ReadIndex - 1);
            Integer nextReadBaseIndex = lookupFromReadIndex.get(coreEndBase.ReadIndex + 1);

            if(prevRefBaseIndex == null || nextRefBaseIndex == null || prevReadBaseIndex == null || nextReadBaseIndex == null)
                return null;

            boolean expandLeft = coreStartBase.RefBase == alignedBases.get(prevRefBaseIndex).RefBase
                    || coreStartBase.ReadBase == alignedBases.get(prevReadBaseIndex).ReadBase;

            boolean expandRight = coreEndBase.RefBase == alignedBases.get(nextRefBaseIndex).RefBase
                    || coreEndBase.ReadBase == alignedBases.get(nextReadBaseIndex).ReadBase;

            if(!expandLeft && !expandRight)
                return new IntPair(coreStart, coreEnd);

            if(expandLeft)
                coreStart--;

            if(expandRight)
                coreEnd++;
        }
    }

    public static int addPadding(final List<AlignedBase> alignedBases, int startIndex, boolean padLeft)
    {
        int inc = padLeft ? -1 : 1;
        for(int i = startIndex + inc; i >= 0 && i < alignedBases.size(); i += inc)
        {
            if(alignedBases.get(i).isExactMatch())
                return i;
        }

        return INVALID_INDEX;
    }

    public static int findFlankBoundary(final List<AlignedBase> alignedBases, int coreBoundaryIndex, boolean leftFlank, int flankSize, boolean inAppendMode)
    {
        int inc = leftFlank ? -1 : 1;
        int readBasesSeen = 0;
        int flankBoundary;
        for(flankBoundary = coreBoundaryIndex + inc; flankBoundary >= 0 && flankBoundary < alignedBases.size(); flankBoundary += inc)
        {
            if(alignedBases.get(flankBoundary).ReadBase != MISSING_BASE)
                readBasesSeen++;

            if(readBasesSeen == flankSize)
                break;
        }

        if(inAppendMode && (flankBoundary < 0 || flankBoundary >= alignedBases.size()))
            flankBoundary = Math.min(Math.max(0, flankBoundary), alignedBases.size() - 1);
        else if(flankBoundary < 0 || flankBoundary >= alignedBases.size())
            return INVALID_INDEX;

        // push alignment out if the flank boundary is in an indel
        while(flankBoundary >= 0 && flankBoundary < alignedBases.size() && alignedBases.get(flankBoundary).isIndel())
            flankBoundary += inc;

        if(inAppendMode && (flankBoundary < 0 || flankBoundary >= alignedBases.size()))
            flankBoundary = Math.min(Math.max(0, flankBoundary), alignedBases.size() - 1);
        else if(flankBoundary < 0 || flankBoundary >= alignedBases.size())
            return INVALID_INDEX;

        return flankBoundary;
    }
}
