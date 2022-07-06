package com.hartwig.hmftools.sage.common;

import static java.lang.Math.max;
import static java.lang.Math.min;

import static com.hartwig.hmftools.sage.SageConstants.CORE_LOW_QUAL_MISMATCH_BASE_LENGTH;
import static com.hartwig.hmftools.sage.SageConstants.MATCHING_BASE_QUALITY;
import static com.hartwig.hmftools.sage.common.ReadContextMatch.CORE;
import static com.hartwig.hmftools.sage.common.ReadContextMatch.FULL;
import static com.hartwig.hmftools.sage.common.ReadContextMatch.NONE;
import static com.hartwig.hmftools.sage.common.ReadContextMatch.PARTIAL;

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

    public String coreString()
    {
        return Bases.length == 0 ? Strings.EMPTY : new String(Bases, LeftCoreIndex, coreLength());
    }

    public String fullString() { return Bases.length == 0 ? Strings.EMPTY : new String(Bases, LeftFlankIndex, length()); }

    public String leftFlankString()
    {
        return Bases.length == 0 ? Strings.EMPTY : new String(Bases, LeftFlankIndex, leftFlankLength());
    }

    public String rightFlankString()
    {
        int rightFlankLength = rightFlankLength();
        return rightFlankLength == 0 ? Strings.EMPTY : new String(Bases, RightCoreIndex + 1, rightFlankLength);
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
    public int coreLength()
    {
        return RightCoreIndex - LeftCoreIndex + 1;
    }
    private int leftFlankLength() { return LeftCoreIndex - LeftFlankIndex; }
    private int rightFlankLength()
    {
        return RightFlankIndex - RightCoreIndex;
    }

    public int corePositionStart() { return Position - (Index - LeftCoreIndex); }
    public int corePositionEnd() { return Position + RightCoreIndex - Index; }

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

    public boolean isCoreCovered(int otherReadIndex, int otherBasesLength)
    {
        int otherLeftCentreIndex = otherLeftCoreIndex(otherReadIndex);

        if(otherLeftCentreIndex < 0)
            return false;

        int otherRightCentreIndex = otherRightCoreIndex(otherReadIndex);
        return otherRightCentreIndex < otherBasesLength;
    }

    public ReadContextMatch matchAtPosition(final IndexedBases other)
    {
        return getMatchType(other, other.length(),null, false, 0);
    }

    public ReadContextMatch matchAtPosition(
            final IndexedBases other, final byte[] otherBaseQuals, boolean wildcardsInCore, int maxCoreMismatches)
    {
        return getMatchType(other, length(), otherBaseQuals, wildcardsInCore, maxCoreMismatches);
    }

    private ReadContextMatch getMatchType(
            final IndexedBases other, int otherLength, final byte[] otherBaseQuals, boolean wildcardsInCore, int maxCoreMismatches)
    {
        int otherReadIndex = other.Index;

        if(otherReadIndex < 0)
            return NONE;

        final byte[] otherBases = other.Bases;

        boolean centreMatch = coreMatch(otherReadIndex, otherBases, otherBaseQuals, wildcardsInCore, maxCoreMismatches);
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

    protected boolean coreMatch(
            int otherRefIndex, final byte[] otherBases, final byte[] otherBaseQuals, boolean wildcardAllowed, int maxCoreMismatches)
    {
        int otherLeftCentreIndex = otherLeftCoreIndex(otherRefIndex);

        if(otherLeftCentreIndex < 0)
            return false;

        int otherRightCentreIndex = otherRightCoreIndex(otherRefIndex);

        if(otherRightCentreIndex >= otherBases.length)
            return false;

        int permittedLowQualMismatches;

        if(maxCoreMismatches > 0)
        {
            permittedLowQualMismatches = 1;
            --maxCoreMismatches;
        }
        else
        {
            permittedLowQualMismatches = 0;
        }

        for(int i = 0; i < coreLength(); i++)
        {
            if(i > 0 && (i % CORE_LOW_QUAL_MISMATCH_BASE_LENGTH) == 0)
            {
                // allow at most 1 mis-match per X bases, so reset it if required
                if(maxCoreMismatches > 0)
                {
                    permittedLowQualMismatches = 1;
                    --maxCoreMismatches;
                }
                else
                {
                    permittedLowQualMismatches = 0;
                }
            }

            int readIndex = LeftCoreIndex + i;
            int otherByteIndex = otherLeftCentreIndex + i;
            byte otherByte = otherBases[otherByteIndex];

            if(Bases[readIndex] == otherByte)
                continue;

            if(wildcardAllowed && otherByte == MATCH_WILDCARD)
                continue;

            if(permittedLowQualMismatches > 0 && isLowBaseQual(otherBaseQuals, otherByteIndex))
            {
                --permittedLowQualMismatches;
                continue;
            }

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

    public String toString()
    {
        if(Bases.length == 0)
            return "";

        return String.format("%s indices(%d-%d-%d-%d-%d)",
                fullString(), LeftFlankIndex, LeftCoreIndex, Index, RightCoreIndex, RightFlankIndex);
    }

}
