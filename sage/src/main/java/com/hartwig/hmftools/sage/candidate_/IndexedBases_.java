package com.hartwig.hmftools.sage.candidate_;

import static java.lang.Math.max;
import static java.lang.Math.min;
import static java.lang.String.format;

import static com.hartwig.hmftools.sage.SageConstants.CORE_LOW_QUAL_MISMATCH_BASE_LENGTH;
import static com.hartwig.hmftools.sage.SageConstants.MATCHING_BASE_QUALITY;
import static com.hartwig.hmftools.sage.common.ReadContextMatch.CORE;
import static com.hartwig.hmftools.sage.common.ReadContextMatch.CORE_PARTIAL;
import static com.hartwig.hmftools.sage.common.ReadContextMatch.FULL;
import static com.hartwig.hmftools.sage.common.ReadContextMatch.NONE;
import static com.hartwig.hmftools.sage.common.ReadContextMatch.PARTIAL;

import com.hartwig.hmftools.sage.common.ReadContextMatch;
import com.hartwig.hmftools.sage.evidence.ReadIndexBases;

import org.apache.commons.lang3.NotImplementedException;
import org.apache.logging.log4j.util.Strings;

/**
 * Represents the bases an indices of a ReadContext.
 * <p>
 * Index - is the index of bases corresponding to position.
 * LeftFlankIndex - this is clipped at 0.
 * RightFlankIndex - this is clipped at the end of bases.
 */
public class IndexedBases_
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

    public IndexedBases_(final int position, final int index, final byte[] bases)
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

    public IndexedBases_(
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

    /**
     * Flanks + core.
     */
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

    /**
     * Index relative to core start.
     */
    public int indexInCore()
    {
        return Index - LeftCoreIndex;
    }

    /**
     * Genome position to bases array index.
     */
    public int index(int position)
    {
        return position - Position + Index;
    }

    public boolean containsPosition(int position)
    {
        int index = index(position);
        return index >= 0 && index < Bases.length;
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

    /**
     * NOTE: Doesn't take into account INDEL bases.
     */
    public int corePositionEnd() { return Position + RightCoreIndex - Index; }

    public byte base(int position) { return Bases[index(position)]; }

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

    /**
     * LeftCoreIndex relative to otherRefIndex.
     */
    private int otherLeftCoreIndex(int otherRefIndex) { return otherRefIndex - Index + LeftCoreIndex; }

    /**
     * RightCoreIndex relative to otherRefIndex.
     */
    private int otherRightCoreIndex(int otherRefIndex)
    {
        return otherRefIndex - Index + RightCoreIndex;
    }

    /**
     * Maps the index to another read index and returns whether the mapped left and right core indices are within [0, otherBasesLength).
     */
    public boolean isCoreCovered(int otherReadIndex, int otherBasesLength)
    {
        // Get left and right core index relative to other read.
        int otherLeftCentreIndex = otherLeftCoreIndex(otherReadIndex);
        int otherRightCentreIndex = otherRightCoreIndex(otherReadIndex);

        return otherLeftCentreIndex >= 0 && otherRightCentreIndex < otherBasesLength;
    }

    /**
     * What type of match does IndexedBases_ have against this.
     * <p>
     * No mismatches allowed, and no wildcard matches allowed.
     * <p>
     * For extra context see the match methods below.
     */
    public ReadContextMatch matchAtPosition(final IndexedBases_ other)
    {
        return getMatchType(other.Index, other.Bases, other.length(), null, false, 0);
    }

    /**
     * What type of match does ReadIndexBases have against this.
     * <p>
     * Uses length (full flank + core length) based of this.
     */
    public ReadContextMatch matchAtPosition(
            final ReadIndexBases readIndexBases, final byte[] readBaseQuals, boolean wildcardsInCore, int maxCoreMismatches)
    {
        return getMatchType(readIndexBases.Index, readIndexBases.Bases, length(), readBaseQuals, wildcardsInCore, maxCoreMismatches);
    }

    /**
     * Match type with another set of bases.
     * <p>
     * For more context, see the following match methods.
     */
    public ReadContextMatch getMatchType(
            final int otherReadIndex, final byte[] otherBases, int otherLength, final byte[] otherBaseQuals,
            boolean wildcardsInCore, int maxCoreMismatches)
    {
        if(otherReadIndex < 0)
            return NONE;

        ReadContextMatch centreMatch = coreMatch(otherReadIndex, otherBases, otherBaseQuals, wildcardsInCore, maxCoreMismatches);
        if(centreMatch == NONE || centreMatch == CORE_PARTIAL)
            return centreMatch;

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

    /**
     * Does the otherBases match this core.
     * @param otherRefIndex
     * @param otherBases
     * @param otherBaseQuals
     * @param wildcardAllowed The otherBases can contain a MATCH_WILDCARD. Note that this does not increment matchCount.
     * @param maxCoreMismatches The number of core low quality mismatches allowed. Note that only one mismatch is allowed per CORE_LOW_QUAL_MISMATCH_BASE_LENGTH bases.
     * @return Returns NONE it a mismatch is found (ex. maxCoreMismatches). If bases were matched and the remapped core goes out of bounds
     * of otherBases then return CORE_PARTIAL, otherwise return CORE.
     */
    public ReadContextMatch coreMatch(
            int otherRefIndex, final byte[] otherBases, final byte[] otherBaseQuals, boolean wildcardAllowed, int maxCoreMismatches)
    {
        // Convert core indices to otherBases space.
        int otherLeftCentreIndex = otherLeftCoreIndex(otherRefIndex);
        int otherRightCentreIndex = otherRightCoreIndex(otherRefIndex);

        // exit if there's no overlap
        if(otherRightCentreIndex < 0 || otherLeftCentreIndex >= otherBases.length)
            return NONE;

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

        boolean partialCoreOverlap = false;
        int matchedBases = 0;

        for(int i = 0; i < coreLength(); i++)
        {
            if(i > 0 && (i % CORE_LOW_QUAL_MISMATCH_BASE_LENGTH) == 0)
            {
                // allow at most 1 mismatch per X bases, so reset it if required
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

            // We are going over the edge of otherBases.
            if(otherByteIndex < 0 || otherByteIndex >= otherBases.length)
            {
                partialCoreOverlap = true;

                if(otherByteIndex >= otherBases.length)
                    break;
                else
                    continue;
            }

            byte otherByte = otherBases[otherByteIndex];

            if(Bases[readIndex] == otherByte)
            {
                ++matchedBases;
                continue;
            }

            if(wildcardAllowed && otherByte == MATCH_WILDCARD)
                continue;

            if(permittedLowQualMismatches > 0 && isLowBaseQual(otherBaseQuals, otherByteIndex))
            {
                --permittedLowQualMismatches;
                continue;
            }

            return NONE;
        }

        return partialCoreOverlap && matchedBases > 0 ? CORE_PARTIAL : CORE;
    }

    /**
     * Returns whether baseQualities at bqIndex is valid and is less that MATCHING_BASE_QUALITY.
     */
    private static boolean isLowBaseQual(final byte[] baseQualities, int bqIndex)
    {
        if(baseQualities == null)
            return false;

        if(bqIndex < 0 || bqIndex >= baseQualities.length)
            return false;

        return baseQualities[bqIndex] < MATCHING_BASE_QUALITY;
    }

    /**
     * How much of the right flank match for otherBases.
     * <p>
     * Allows wildcard matches.
     *
     * @return -1 or the max of the possible right flank length (note that is may be clipped because of the size of otherBases and the
     * position of otherRefIndex in this).
     */
    public int rightFlankMatchingBases(int otherRefIndex, final byte[] otherBases)
    {
        return rightFlankMatchingBases(otherRefIndex, otherBases, null);
    }

    /**
     * How much of the right flank match for otherBases.
     * <p>
     * Allows low base qual mismatches and wildcard matches.
     *
     * @return -1 or the max possible right flank length (note that is may be clipped because of the size of otherBases and the
     * position of otherRefIndex in this).
     */
    public int rightFlankMatchingBases(int otherRefIndex, byte[] otherBases, final byte[] otherBaseQuals)
    {
        // Maps right centre index using otherRefIndex.
        int otherRightCentreIndex = otherRefIndex + RightCoreIndex - Index;

        // Computers right flank length, clips on the end of otherBases.
        int otherRightFlankLength = min(otherBases.length - 1 - otherRightCentreIndex, FlankSize);
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

    /**
     * How much of the left flank match for otherBases.
     * <p>
     * Allows wildcard mismatches.
     *
     * @return -1 or the max possible left flank length (note that is may be clipped because of the size of otherBases and the
     * position of otherRefIndex in this).
     */
    public int leftFlankMatchingBases(int otherRefIndex, final byte[] otherBases)
    {
        return leftFlankMatchingBases(otherRefIndex, otherBases, null);
    }

    /**
     * How much of the left flank match for otherBases.
     * <p>
     * Allows low base qual mismatches and wildcard mismatches.
     *
     * @return -1 or the max possible left flank length (note that is may be clipped because of the size of otherBases and the
     * position of otherRefIndex in this).
     */
    public int leftFlankMatchingBases(int otherRefIndex, byte[] otherBases, final byte[] otherBaseQuals)
    {
        int otherLeftCentreIndex = otherRefIndex + LeftCoreIndex - Index;
        int otherLeftFlankLength =  min(otherLeftCentreIndex, FlankSize);
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

    @Override
    public String toString()
    {
        if(Bases.length == 0)
            return "";

        return format("indices(%d-%d-%d-%d-%d) corePos(%d - %d - %d) bases(%s - %s - %s)",
                LeftFlankIndex, LeftCoreIndex, Index, RightCoreIndex, RightFlankIndex,
                corePositionStart(), Position, corePositionEnd(), leftFlankString(), coreString(), rightFlankString());
    }
}
