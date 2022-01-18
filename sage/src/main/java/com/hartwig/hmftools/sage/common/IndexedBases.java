package com.hartwig.hmftools.sage.common;

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
        LeftFlankIndex = Math.max(0, leftCoreIndex - flankSize);
        RightFlankIndex = Math.min(bases.length - 1, rightCoreIndex + flankSize);
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
        return Math.min(LeftCoreIndex, Bases.length - RightCoreIndex - 1);
    }

    public byte[] trinucleotideContext(int position)
    {
        return new byte[] { base(position - 1), base(position), base(position + 1) };
    }

    public boolean flanksComplete()
    {
        return leftFlankLength() == FlankSize && rightFlankLength() == FlankSize;
    }

    private int otherLeftCentreIndex(int otherRefIndex) { return otherRefIndex + LeftCoreIndex - Index; }
    private int otherRightCentreIndex(int otherRefIndex)
    {
        return otherRefIndex + RightCoreIndex - Index;
    }

    public boolean isCentreCovered(int otherReadIndex, byte[] otherBases)
    {
        int otherLeftCentreIndex = otherLeftCentreIndex(otherReadIndex);

        if(otherLeftCentreIndex < 0)
            return false;

        int otherRightCentreIndex = otherRightCentreIndex(otherReadIndex);
        return otherRightCentreIndex < otherBases.length;
    }

    public ReadContextMatch matchAtPosition(final IndexedBases other)
    {
        return getMatchType(false, other, null);
    }

    public ReadContextMatch matchAtPosition(final IndexedBases other, boolean wildcardAllowedInCore, final int[] otherBaseQuals)
    {
        return getMatchType(wildcardAllowedInCore, other, otherBaseQuals);
    }

    private ReadContextMatch getMatchType(
            boolean wildcardAllowedInCore, final IndexedBases other, final int[] otherBaseQuals)
    {
        int otherReadIndex = other.Index;

        if(otherReadIndex < 0)
            return NONE;

        final byte[] otherBases = other.Bases;

        int otherCoreBaseQualStartIndex = other.LeftCoreIndex - other.LeftFlankIndex;
        boolean centreMatch = coreMatch(wildcardAllowedInCore, otherReadIndex, otherBases, otherBaseQuals, otherCoreBaseQualStartIndex);
        if(!centreMatch)
            return NONE;

        int otherLeftFlankOffset = other.LeftFlankIndex;
        int leftFlankingBases = leftFlankMatchingBases(otherReadIndex, otherBases, otherBaseQuals, otherLeftFlankOffset);
        if(leftFlankingBases < 0)
            return CORE;

        int rightFlankingBases = rightFlankMatchingBases(otherReadIndex, otherBases, otherBaseQuals, otherLeftFlankOffset);
        if(rightFlankingBases < 0)
            return CORE;

        int leftFlankLength = leftFlankLength();
        int rightFlankLength = rightFlankLength();

        if(leftFlankingBases != leftFlankLength && rightFlankingBases != rightFlankLength)
            return CORE;

        return length() == other.length() && leftFlankingBases == leftFlankLength && rightFlankingBases == rightFlankLength ? FULL : PARTIAL;
    }

    protected boolean coreMatch(boolean wildcardAllowed, int otherRefIndex, final byte[] otherBases)
    {
        return coreMatch(wildcardAllowed, otherRefIndex, otherBases, null, 0);
    }

    protected boolean coreMatch(
            boolean wildcardAllowed, int otherRefIndex, final byte[] otherBases, final int[] otherBaseQuals, int bqStartIndex)
    {
        int otherLeftCentreIndex = otherLeftCentreIndex(otherRefIndex);

        if(otherLeftCentreIndex < 0)
            return false;

        int otherRightCentreIndex = otherRightCentreIndex(otherRefIndex);

        if(otherRightCentreIndex >= otherBases.length)
            return false;

        for(int i = 0; i < coreLength(); i++)
        {
            byte ourByte = Bases[LeftCoreIndex + i];
            byte otherByte = otherBases[otherLeftCentreIndex + i];

            if(!bytesMatch(wildcardAllowed, ourByte, otherByte))
            {
                if(isLowBaseQual(otherBaseQuals, bqStartIndex + i))
                    continue;

                return false;
            }
        }

        return true;
    }

    private static boolean isLowBaseQual(final int[] baseQualities, int bqIndex)
    {
        if(baseQualities == null)
            return false;

        if(bqIndex < 0 || bqIndex >= baseQualities.length)
            return false;

        return baseQualities[bqIndex] < MATCHING_BASE_QUALITY;
    }

    private boolean bytesMatch(final boolean wildcardAllowed, byte ours, byte other)
    {
        return (wildcardAllowed && other == MATCH_WILDCARD) || ours == other;
    }

    protected int rightFlankMatchingBases(int otherRefIndex, final byte[] otherBases)
    {
        return rightFlankMatchingBases(otherRefIndex, otherBases, null, 0);
    }

    protected int rightFlankMatchingBases(int otherRefIndex, byte[] otherBases, final int[] otherBaseQuals, int bqLeftFlankOffset)
    {
        int otherRightCentreIndex = otherRefIndex + RightCoreIndex - Index;
        int otherRightFlankLength = Math.min(otherBases.length - 1, otherRightCentreIndex + FlankSize) - otherRightCentreIndex;
        int maxLength = Math.min(rightFlankLength(), otherRightFlankLength);

        for(int i = 1; i <= maxLength; i++)
        {
            int otherByteIndex = otherRightCentreIndex + i;
            byte otherByte = otherBases[otherByteIndex];

            if(Bases[RightCoreIndex + i] != otherByte && otherByte != MATCH_WILDCARD)
            {
                if(isLowBaseQual(otherBaseQuals, otherByteIndex - bqLeftFlankOffset))
                    continue;

                return -1;
            }
        }

        return maxLength;
    }

    protected int leftFlankMatchingBases(int otherRefIndex, final byte[] otherBases)
    {
        return leftFlankMatchingBases(otherRefIndex, otherBases, null, 0);
    }

    protected int leftFlankMatchingBases(int otherRefIndex, byte[] otherBases, final int[] otherBaseQuals, int bqLeftFlankOffset)
    {
        int otherLeftCentreIndex = otherRefIndex + LeftCoreIndex - Index;
        int otherLeftFlankLength = otherLeftCentreIndex - Math.max(0, otherLeftCentreIndex - FlankSize);
        int totalLength = Math.min(leftFlankLength(), otherLeftFlankLength);

        for(int i = 1; i <= totalLength; i++)
        {
            int otherByteIndex = otherLeftCentreIndex - i;
            byte otherByte = otherBases[otherLeftCentreIndex - i];

            if(Bases[LeftCoreIndex - i] != otherByte && otherByte != MATCH_WILDCARD)
            {
                if(isLowBaseQual(otherBaseQuals, otherByteIndex - bqLeftFlankOffset))
                    continue;

                return -1;
            }
        }

        return totalLength;
    }

    public boolean phased(int offset, final IndexedBases other)
    {
        int otherReadIndex = other.Index + offset;

        boolean centreMatch = coreMatch(false, otherReadIndex, other.Bases, null, 0);

        if(!centreMatch)
            return false;

        boolean otherCentreMatch = other.coreMatch(false, Index - offset, Bases, null, 0);

        if(!otherCentreMatch)
            return false;

        int leftFlankingBases = leftFlankMatchingBases(otherReadIndex, other.Bases, null, 0);

        if(leftFlankingBases < 0)
            return false;

        int rightFlankingBases = rightFlankMatchingBases(otherReadIndex, other.Bases, null, 0);

        return rightFlankingBases >= 0 && (rightFlankingBases >= FlankSize || leftFlankingBases >= FlankSize);
    }

    public static IndexedBases resize(
            final int position, final int recordIndex, final int recordLeftCoreIndex,
            final int recordRightCoreIndex, final int flankSize, final int additionalFlank, final byte[] recordBases)
    {
        int recordLeftFlankIndex = Math.max(0, recordLeftCoreIndex - flankSize - additionalFlank);
        int recordLeftFlankLength = recordLeftCoreIndex - recordLeftFlankIndex;
        int recordRightFlankIndex = Math.min(recordBases.length - 1, recordRightCoreIndex + flankSize + additionalFlank);

        int rightCentreIndex = recordLeftFlankLength + recordRightCoreIndex - recordLeftCoreIndex;
        int index = recordLeftFlankLength + recordIndex - recordLeftCoreIndex;
        byte[] bases = Arrays.copyOfRange(recordBases, recordLeftFlankIndex, recordRightFlankIndex + 1);
        return new IndexedBases(position, index, recordLeftFlankLength, rightCentreIndex, flankSize, bases);
    }

}
