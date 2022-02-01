package com.hartwig.hmftools.sage.common;

import static java.lang.Math.max;
import static java.lang.Math.min;

import static com.hartwig.hmftools.common.utils.sv.BaseRegion.positionWithin;
import static com.hartwig.hmftools.common.utils.sv.BaseRegion.positionsOverlap;
import static com.hartwig.hmftools.sage.SageConstants.MATCHING_BASE_QUALITY;
import static com.hartwig.hmftools.sage.common.ReadContextMatch.CORE;
import static com.hartwig.hmftools.sage.common.ReadContextMatch.FULL;
import static com.hartwig.hmftools.sage.common.ReadContextMatch.NONE;
import static com.hartwig.hmftools.sage.common.ReadContextMatch.PARTIAL;

import java.util.Arrays;

import org.apache.logging.log4j.util.Strings;

public class IndexedBases
{
    public final int Position;
    public final int Index;
    public final int FlankSize;
    public final int LeftFlankIndex;
    public final int LeftCoreIndex;
    public final int RightCoreIndex;
    public final int RightFlankIndex;
    public final byte[] Bases;

    public static final byte MATCH_WILDCARD = (byte) '.';

    public IndexedBases(final int position, final int index, final byte[] bases)
    {
        Position = position;
        Index = index;
        LeftCoreIndex = index;
        RightCoreIndex = index;
        Bases = bases;
        LeftFlankIndex = index;
        RightFlankIndex = index;
        FlankSize = 0;
    }

    public IndexedBases(
            final int position, final int index, final int leftCoreIndex, int rightCoreIndex, int flankSize, final byte[] bases)
    {
        Position = position;
        Index = index;
        LeftCoreIndex = leftCoreIndex;
        RightCoreIndex = rightCoreIndex;
        Bases = bases;
        LeftFlankIndex = max(0, leftCoreIndex - flankSize);
        RightFlankIndex = min(bases.length - 1, rightCoreIndex + flankSize);
        FlankSize = flankSize;
    }

    public String centerString()
    {
        return Bases.length == 0 ? Strings.EMPTY : new String(Bases, LeftCoreIndex, coreLength());
    }

    public String leftFlankString()
    {
        return Bases.length == 0 ? Strings.EMPTY : new String(Bases, LeftFlankIndex, leftFlankLength());
    }

    public String rightFlankString()
    {
        int rightFlankLength = rightFlankLength();
        return rightFlankLength == 0 ? Strings.EMPTY : new String(Bases, RightCoreIndex + 1, rightFlankLength);
    }

    @Override
    public String toString()
    {
        return Bases.length == 0 ? Strings.EMPTY : new String(Bases, LeftFlankIndex, length());
    }

    public int indexInCore()
    {
        return Index - LeftCoreIndex;
    }

    public int index(int position)
    {
        return position - Position + Index;
    }

    public int length() { return RightFlankIndex - LeftFlankIndex + 1; }
    private int coreLength()
    {
        return RightCoreIndex - LeftCoreIndex + 1;
    }
    private int leftFlankLength() { return LeftCoreIndex - LeftFlankIndex; }
    private int rightFlankLength()
    {
        return RightFlankIndex - RightCoreIndex;
    }

    public byte base(int position)
    {
        return Bases[position - Position + Index];
    }

    public int maxFlankLength()
    {
        return min(LeftCoreIndex, Bases.length - RightCoreIndex - 1);
    }

    public byte[] trinucleotideContext(int position)
    {
        return new byte[] { base(position - 1), base(position), base(position + 1) };
    }

    public boolean flanksComplete()
    {
        return leftFlankLength() == FlankSize && rightFlankLength() == FlankSize;
    }

    // calculate implied coords for another IndexedBase object
    private int otherLeftCoreIndex(int otherRefIndex) { return otherRefIndex - Index + LeftCoreIndex; }
    private int otherRightCoreIndex(int otherRefIndex)
    {
        return otherRefIndex- Index + RightCoreIndex;
    }

    public boolean isCentreCovered(int otherReadIndex, byte[] otherBases)
    {
        int otherLeftCentreIndex = otherLeftCoreIndex(otherReadIndex);

        if(otherLeftCentreIndex < 0)
            return false;

        int otherRightCentreIndex = otherRightCoreIndex(otherReadIndex);
        return otherRightCentreIndex < otherBases.length;
    }

    public ReadContextMatch matchAtPosition(final IndexedBases other)
    {
        return getMatchType(false, other, other.length(),null, 0);
    }

    public ReadContextMatch matchAtPosition(
            final IndexedBases other, boolean wildcardAllowedInCore, final byte[] otherBaseQuals, int nonIndelLength)
    {
        return getMatchType(wildcardAllowedInCore, other, length(), otherBaseQuals, nonIndelLength);
    }

    private ReadContextMatch getMatchType(
            boolean wildcardAllowedInCore, final IndexedBases other, int otherLength, final byte[] otherBaseQuals, int nonIndelLength)
    {
        int otherReadIndex = other.Index;

        if(otherReadIndex < 0)
            return NONE;

        final byte[] otherBases = other.Bases;

        boolean centreMatch = coreMatch(wildcardAllowedInCore, otherReadIndex, otherBases);
        if(!centreMatch)
            return NONE;

        int leftFlankingBases = leftFlankMatchingBases(otherReadIndex, otherBases, otherBaseQuals);
        if(leftFlankingBases < 0)
            return CORE;

        int rightFlankingBases = rightFlankMatchingBases(otherReadIndex, otherBases, otherBaseQuals);
        if(rightFlankingBases < 0)
            return CORE;

        int leftFlankLength = leftFlankLength();
        int rightFlankLength = rightFlankLength();

        if(leftFlankingBases != leftFlankLength && rightFlankingBases != rightFlankLength)
            return CORE;

        return length() == otherLength && leftFlankingBases == leftFlankLength && rightFlankingBases == rightFlankLength ? FULL : PARTIAL;
    }

    protected boolean coreMatch(boolean wildcardAllowed, int otherRefIndex, final byte[] otherBases)
    {
        int otherLeftCentreIndex = otherLeftCoreIndex(otherRefIndex);

        if(otherLeftCentreIndex < 0)
            return false;

        int otherRightCentreIndex = otherRightCoreIndex(otherRefIndex);

        if(otherRightCentreIndex >= otherBases.length)
            return false;

        for(int i = 0; i < coreLength(); i++)
        {
            int readIndex = LeftCoreIndex + i;
            byte ourByte = Bases[readIndex];
            byte otherByte = otherBases[otherLeftCentreIndex + i];

            if(!bytesMatch(wildcardAllowed, ourByte, otherByte))
                return false;
        }

        return true;
    }

    private static boolean isLowBaseQual(final byte[] baseQualities, int bqIndex)
    {
        if(baseQualities == null)
            return false;

        if(bqIndex < 0 || bqIndex >= baseQualities.length)
            return false;

        return baseQualities[bqIndex] < MATCHING_BASE_QUALITY;
    }

    private static boolean bytesMatch(final boolean wildcardAllowed, byte ours, byte other)
    {
        return (wildcardAllowed && other == MATCH_WILDCARD) || ours == other;
    }

    protected int rightFlankMatchingBases(int otherRefIndex, final byte[] otherBases)
    {
        return rightFlankMatchingBases(otherRefIndex, otherBases, null);
    }

    protected int rightFlankMatchingBases(int otherRefIndex, byte[] otherBases, final byte[] otherBaseQuals)
    {
        int otherRightCentreIndex = otherRefIndex + RightCoreIndex - Index;
        int otherRightFlankLength = min(otherBases.length - 1, otherRightCentreIndex + FlankSize) - otherRightCentreIndex;
        int maxLength = min(rightFlankLength(), otherRightFlankLength);

        for(int i = 1; i <= maxLength; i++)
        {
            int otherByteIndex = otherRightCentreIndex + i;
            byte otherByte = otherBases[otherByteIndex];

            if(Bases[RightCoreIndex + i] != otherByte && otherByte != MATCH_WILDCARD)
            {
                if(isLowBaseQual(otherBaseQuals, otherByteIndex))
                    continue;

                return -1;
            }
        }

        return maxLength;
    }

    protected int leftFlankMatchingBases(int otherRefIndex, final byte[] otherBases)
    {
        return leftFlankMatchingBases(otherRefIndex, otherBases, null);
    }

    protected int leftFlankMatchingBases(int otherRefIndex, byte[] otherBases, final byte[] otherBaseQuals)
    {
        int otherLeftCentreIndex = otherRefIndex + LeftCoreIndex - Index;
        int otherLeftFlankLength = otherLeftCentreIndex - max(0, otherLeftCentreIndex - FlankSize);
        int totalLength = min(leftFlankLength(), otherLeftFlankLength);

        for(int i = 1; i <= totalLength; i++)
        {
            int otherByteIndex = otherLeftCentreIndex - i;
            byte otherByte = otherBases[otherLeftCentreIndex - i];

            if(Bases[LeftCoreIndex - i] != otherByte && otherByte != MATCH_WILDCARD)
            {
                if(isLowBaseQual(otherBaseQuals, otherByteIndex))
                    continue;

                return -1;
            }
        }

        return totalLength;
    }

    public boolean phasedNew(int offset, final IndexedBases other)
    {
        int otherRefIndex = other.Index + offset;

        int otherLeftCoreIndex = otherLeftCoreIndex(otherRefIndex);

        if(otherLeftCoreIndex < 0)
            return false;

        int otherRightCoreIndex = otherRightCoreIndex(otherRefIndex);

        if(otherRightCoreIndex >= other.Bases.length)
            return false;

        // test this core vs the other variant's bases
        int coreLength = RightCoreIndex - LeftCoreIndex + 1;

        boolean hasCoreMatch = true;

        for(int i = 0; i < coreLength; i++)
        {
            int readIndex = LeftCoreIndex + i;
            byte ourByte = Bases[readIndex];

            int otherReadIndex = otherLeftCoreIndex + i;
            byte otherByte = other.Bases[otherReadIndex];

            if(ourByte != otherByte)
            {
                if(positionWithin(otherReadIndex, other.LeftFlankIndex, other.RightFlankIndex))
                    return false;

                hasCoreMatch = false;
            }
        }

        if(!hasCoreMatch)
        {
            // test other core vs this variant's bases
            hasCoreMatch = true;

            int adjRefIndex = Index - offset;
            int adjLeftCoreIndex = other.otherLeftCoreIndex(adjRefIndex);

            if(adjLeftCoreIndex < 0)
                return false;

            int adjRightCoreIndex = other.otherRightCoreIndex(adjRefIndex);

            if(adjRightCoreIndex >= Bases.length)
                return false;

            for(int i = 0; i < other.coreLength(); i++)
            {
                int otherReadIndex = other.LeftCoreIndex + i;
                byte otherByte = other.Bases[otherReadIndex];

                int readIndex = adjLeftCoreIndex + i;
                byte ourByte = Bases[readIndex];

                if(ourByte != otherByte)
                {
                    if(positionWithin(readIndex, LeftFlankIndex, RightFlankIndex))
                        return false;

                    hasCoreMatch = false;
                }
            }
        }

        if(!hasCoreMatch) // neither variant has a core match in their other's read bases
            return false;

        // test left flank bases

        // int otherLeftCentreIndex = otherRefIndex + LeftCoreIndex - Index;
        int otherLeftFlankLength = otherLeftCoreIndex - max(0, otherLeftCoreIndex - FlankSize);
        int leftFlankLength = LeftCoreIndex - LeftFlankIndex;

        int totalLength = min(leftFlankLength, otherLeftFlankLength);

        for(int i = 1; i <= totalLength; i++)
        {
            byte ourByte = Bases[LeftCoreIndex - i];

            int otherReadIndex = otherLeftCoreIndex - i;
            byte otherByte = other.Bases[otherReadIndex];

            if(ourByte != otherByte && positionWithin(otherReadIndex, LeftFlankIndex, RightFlankIndex))
                return false;
        }

        // test right flank bases

        // int otherRightCentreIndex = otherRefIndex + RightCoreIndex - Index;
        int otherRightFlankLength = min(other.Bases.length - 1, otherRightCoreIndex + FlankSize) - otherRightCoreIndex;
        int rightFlankLength = RightFlankIndex - RightCoreIndex;
        int maxLength = min(rightFlankLength, otherRightFlankLength);

        for(int i = 1; i <= maxLength; i++)
        {
            byte ourByte = Bases[RightCoreIndex + i];
            int otherReadIndex = otherRightCoreIndex + i;
            byte otherByte = other.Bases[otherReadIndex];

            if(ourByte != otherByte &&  positionWithin(otherReadIndex, LeftFlankIndex, RightFlankIndex))
                return false;
        }

        return true;
    }

    public boolean phased(int offset, final IndexedBases other)
    {
        int otherReadIndex = other.Index + offset;

        boolean centreMatch = coreMatch(false, otherReadIndex, other.Bases);

        if(!centreMatch)
            return false;

        boolean otherCentreMatch = other.coreMatch(false, Index - offset, Bases);

        if(!otherCentreMatch)
            return false;

        int leftFlankingBases = leftFlankMatchingBases(otherReadIndex, other.Bases);

        if(leftFlankingBases < 0)
            return false;

        int rightFlankingBases = rightFlankMatchingBases(otherReadIndex, other.Bases);

        return rightFlankingBases >= 0 && (rightFlankingBases >= FlankSize || leftFlankingBases >= FlankSize);
    }

    public boolean phasedNewByPos(int offset, final IndexedBases other)
    {
        int otherAdjReadIndex = other.Index + offset;
        int otherAdjPosition = Position;

        int posIndexOffset = Position - Index;
        int coreStart = Position - (Index - LeftCoreIndex);
        int coreEnd = Position + (RightCoreIndex - Index);
        int leftFlankPos = Position - (Index - LeftFlankIndex);
        int rightFlankPos = Position + (RightFlankIndex - Index);

        int otherPosIndexOffset = otherAdjPosition - otherAdjReadIndex;
        int otherCoreStart = otherAdjPosition - (otherAdjReadIndex - other.LeftCoreIndex);
        int otherCoreEnd = otherAdjPosition + (other.RightCoreIndex - otherAdjReadIndex);
        int otherLeftFlankPos = otherAdjPosition - (otherAdjReadIndex - other.LeftFlankIndex);
        int otherRightFlankPos = otherAdjPosition + (other.RightFlankIndex - otherAdjReadIndex);

        /*
        int otherPosIndexOffset = other.Position - other.Index;
        int otherCoreStart = other.Position - (other.Index - other.LeftCoreIndex);
        int otherCoreEnd = other.Position + (other.RightCoreIndex - other.Index);
        int otherLeftFlankPos = other.Position - (other.Index - other.LeftFlankIndex);
        int otherRightFlankPos = other.Position + (other.RightFlankIndex - other.Index);
        */

        // determine overlap of core + flanks - these must exact match
        if(positionsOverlap(leftFlankPos, rightFlankPos, otherLeftFlankPos, otherRightFlankPos))
        {
            int overlapStart = max(leftFlankPos, otherLeftFlankPos);
            int overlapEnd = min(rightFlankPos, otherRightFlankPos);

            for(int pos = overlapStart; pos <= overlapEnd; ++pos)
            {
                int readIndex = pos - posIndexOffset;
                byte thisByte = Bases[readIndex];

                int otherReadIndex = pos - otherPosIndexOffset;
                byte otherByte = other.Bases[otherReadIndex];

                if(thisByte != otherByte)
                    return false;
            }
        }

        // additionally check any part of core that is only overlapped by the extended flanks - this must also match
        if(coreEnd > otherRightFlankPos)
        {
            if(!checkCoreInExtendedBasesRight(coreStart, coreEnd, Bases, posIndexOffset, otherRightFlankPos, other.Bases, otherPosIndexOffset))
                return false;
        }
        else if(otherCoreEnd > rightFlankPos)
        {
            if(!checkCoreInExtendedBasesRight(otherCoreStart, otherCoreEnd, other.Bases, otherPosIndexOffset, rightFlankPos, Bases, posIndexOffset))
                return false;
        }

        if(coreStart < otherLeftFlankPos)
        {
            if(!checkCoreInExtendedBasesLeft(coreStart, coreEnd,  Bases, posIndexOffset, otherLeftFlankPos, other.Bases, otherPosIndexOffset))
                return false;
        }
        else if(otherCoreStart < leftFlankPos)
        {
            if(!checkCoreInExtendedBasesLeft(otherCoreStart, otherCoreEnd, other.Bases, otherPosIndexOffset, leftFlankPos, Bases, posIndexOffset))
                return false;
        }

        return true;
    }

    private static boolean checkCoreInExtendedBasesRight(
            int coreStart, int coreEnd, final byte[] bases, int posIndexOffset, int otherRightFlankPos, final byte[] otherBases, int otherPosIndexOffset)
    {
        int startPos = max(coreStart, otherRightFlankPos + 1);

        for(int pos = startPos; pos <= coreEnd; ++pos)
        {
            int readIndex = pos - posIndexOffset;
            byte thisByte = bases[readIndex];

            int otherReadIndex = pos - otherPosIndexOffset;

            if(otherReadIndex >= otherBases.length)
                return false;

            byte otherByte = otherBases[otherReadIndex];

            if(thisByte != otherByte)
                return false;
        }

        return true;
    }

    private static boolean checkCoreInExtendedBasesLeft(
            int coreStart, int coreEnd, final byte[] bases, int posIndexOffset, int otherLeftFlankPos, final byte[] otherBases, int otherPosIndexOffset)
    {
        int endPos = min(coreEnd, otherLeftFlankPos - 1);

        for(int pos = coreStart; pos <= endPos; ++pos)
        {
            int readIndex = pos - posIndexOffset;
            byte thisByte = bases[readIndex];

            int otherReadIndex = pos - otherPosIndexOffset;

            if(otherReadIndex < 0)
                return false;

            byte otherByte = otherBases[otherReadIndex];

            if(thisByte != otherByte)
                return false;
        }

        return true;
    }

    public static IndexedBases resize(
            final int position, final int recordIndex, final int recordLeftCoreIndex,
            final int recordRightCoreIndex, final int flankSize, final int additionalFlank, final byte[] recordBases)
    {
        int recordLeftFlankIndex = max(0, recordLeftCoreIndex - flankSize - additionalFlank);
        int recordLeftFlankLength = recordLeftCoreIndex - recordLeftFlankIndex;
        int recordRightFlankIndex = min(recordBases.length - 1, recordRightCoreIndex + flankSize + additionalFlank);

        int rightCentreIndex = recordLeftFlankLength + recordRightCoreIndex - recordLeftCoreIndex;
        int index = recordLeftFlankLength + recordIndex - recordLeftCoreIndex;
        byte[] bases = Arrays.copyOfRange(recordBases, recordLeftFlankIndex, recordRightFlankIndex + 1);
        return new IndexedBases(position, index, recordLeftFlankLength, rightCentreIndex, flankSize, bases);
    }

}
