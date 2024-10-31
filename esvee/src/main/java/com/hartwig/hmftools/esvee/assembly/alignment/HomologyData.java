package com.hartwig.hmftools.esvee.assembly.alignment;

import static java.lang.Math.abs;
import static java.lang.Math.max;
import static java.lang.String.format;

import static com.hartwig.hmftools.esvee.assembly.alignment.BreakendBuilder.segmentOrientation;

import com.hartwig.hmftools.common.codon.Nucleotides;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeInterface;
import com.hartwig.hmftools.common.genome.region.Orientation;

public class HomologyData
{
    public final String Homology;
    public final int ExactStart;
    public final int ExactEnd;
    public final int InexactStart;
    public final int InexactEnd;

    public HomologyData(final String homology, final int exactStart, final int exactEnd, final int inexactStart, final int inexactEnd)
    {
        Homology = homology;
        ExactStart = exactStart;
        ExactEnd = exactEnd;
        InexactStart = inexactStart;
        InexactEnd = inexactEnd;
    }

    public String toString() { return format("%s exact(%d,%d) inexact(%d,%d)",
            Homology.isEmpty() ? "none" : Homology, ExactStart, ExactEnd, InexactStart, InexactEnd); }

    public boolean isSymmetrical() { return abs(InexactStart) == InexactEnd; }

    public int positionAdjustment(final Orientation orientation)
    {
        return orientation.isForward() ? -InexactEnd : abs(InexactStart);
    }

    public HomologyData invert(final boolean reversePositions, final boolean reverseBases)
    {
        return new HomologyData(
                reverseBases ? Nucleotides.reverseComplementBases(Homology) : Homology,
                reversePositions ? -ExactEnd : ExactStart, reversePositions ? abs(ExactStart) : ExactEnd,
                reversePositions ? -InexactEnd : InexactStart, reversePositions ? abs(InexactStart) : InexactEnd);
    }

    public static HomologyData determineHomology(final String fullSequence, final AlignData leftAlignment, final AlignData rightAlignment)
    {
        int overlapLength = leftAlignment.sequenceEnd() - rightAlignment.sequenceStart() + 1;

        if(overlapLength <= 0)
            return null;

        int leftIndexEnd = leftAlignment.segmentLength() - 1;
        int leftIndexStart = leftIndexEnd - overlapLength + 1;

        int rightIndexStart = 0;
        int rightIndexEnd = overlapLength - 1;

        String overlapBases = fullSequence.substring(rightAlignment.sequenceStart(), leftAlignment.sequenceEnd() + 1);

        if(leftAlignment.isReverse())
            overlapBases = Nucleotides.reverseComplementBases(overlapBases);

        MdTag leftMdTag = new MdTag(leftAlignment.mdTag());
        MdTag rightMdTag = new MdTag(rightAlignment.mdTag());

        if(!leftMdTag.hasMismatches() && !rightMdTag.hasMismatches())
        {
            int halfOverlap = overlapLength / 2;
            int exactStart = (overlapLength % 2) == 0 ? halfOverlap : halfOverlap + 1; // round up if an odd length
            int exactEnd = overlapLength - exactStart;
            return new HomologyData(overlapBases, -exactStart, exactEnd, -exactStart, exactEnd);
        }

        byte[] leftMdSeq = leftMdTag.extractSubSequence(leftIndexStart, leftIndexEnd, leftAlignment.isReverse());
        byte[] rightMdSeq = rightMdTag.extractSubSequence(rightIndexStart, rightIndexEnd, rightAlignment.isReverse());

        if(leftMdSeq.length != rightMdSeq.length)
            return null;

        // count mismatches between the 2 sequence arrays progressively along the full sequence
        int[] indexMismatches = new int[leftMdSeq.length + 1];
        int minMismatchCount = -1;
        int minMismatchEndIndex = -1;

        for(int endIndex = 0; endIndex <= leftMdSeq.length; ++endIndex)
        {
            int mismatchCount = 0;

            for(int i = 0; i < leftMdSeq.length; ++i)
            {
                if(i < endIndex)
                {
                    if(leftMdSeq[i] != MdTag.MATCH_BYTE)
                        ++mismatchCount;
                }
                else
                {
                    if(rightMdSeq[i] != MdTag.MATCH_BYTE)
                        ++mismatchCount;
                }
            }

            indexMismatches[endIndex] = mismatchCount;

            if(minMismatchCount == -1 || mismatchCount <= minMismatchCount)
            {
                minMismatchCount = mismatchCount;
                minMismatchEndIndex = endIndex;
            }
        }

        // first a simple check for full exact homology
        if(minMismatchCount == 0 && minMismatchEndIndex == overlapLength - 1)
        {
            int halfOverlap = overlapLength / 2;
            int exactStart = (overlapLength % 2) == 0 ? halfOverlap : halfOverlap + 1; // round up if an odd length
            int exactEnd = overlapLength - exactStart;
            return new HomologyData(overlapBases, -exactStart, exactEnd, -exactStart, exactEnd);
        }

        // find the range with minimum mismatches
        int rangeStart = -1;
        int rangeEnd = -1;
        int longestRangeStart = -1;
        int longestRangeEnd = -1;

        for(int i = 0; i < indexMismatches.length; ++i)
        {
            if(indexMismatches[i] == minMismatchCount)
            {
                if(rangeStart < 0)
                    rangeStart = i;
            }
            else if(rangeStart >= 0)
            {
                rangeEnd = i - 1;

                if(longestRangeStart < 0 || (rangeEnd - rangeStart + 1 > longestRangeEnd - longestRangeStart + 1))
                {
                    longestRangeStart = rangeStart;
                    longestRangeEnd = rangeEnd;
                }

                // reset
                rangeStart = -1;
                rangeEnd = -1;
            }
        }

        if(rangeStart >= 0 && rangeEnd < 0)
        {
            // extended until the end
            longestRangeStart = rangeStart;
            longestRangeEnd = indexMismatches.length - 1;
        }

        int rangeLength = longestRangeEnd - longestRangeStart + 1;
        int halfRange = rangeLength / 2;

        int exactStart = 0;
        int exactEnd = 0;
        String exactHomology = "";

        if(longestRangeStart >= 0 && longestRangeEnd > longestRangeStart)
        {
            exactStart = max(halfRange, 1); // round up if an odd length
            exactEnd = rangeLength - exactStart - 1;

            exactHomology = overlapBases.substring(longestRangeStart, longestRangeEnd);
        }

        int inexactStart = longestRangeStart + halfRange;
        int inexactEnd = overlapLength - inexactStart;

        return new HomologyData(exactHomology, -exactStart, exactEnd, -inexactStart, inexactEnd);
    }

    public boolean matches(final HomologyData other)
    {
        return other.Homology.equals(Homology) && ExactStart == other.ExactStart && ExactEnd == other.ExactEnd
                && InexactStart == other.InexactStart && InexactEnd == other.InexactEnd;
    }

    public static HomologyData determineIndelHomology(final String refBasesStart, final String refBasesEnd, final int overlap)
    {
        if(refBasesStart.length() != refBasesEnd.length())
            return null;

        // first a simple check for full exact homology
        if(refBasesStart.equals(refBasesEnd))
        {
            int halfOverlap = overlap / 2;
            int exactStart = (overlap % 2) == 0 ? halfOverlap : halfOverlap + 1; // round up if an odd length
            int exactEnd = overlap - exactStart;
            return new HomologyData(refBasesStart, -exactStart, exactEnd, -exactStart, exactEnd);
        }

        int exactMatch = 0;
        StringBuilder sb = new StringBuilder();

        for(int i = 0; i < overlap; ++i)
        {
            if(refBasesStart.charAt(i) == refBasesEnd.charAt(i))
            {
                sb.append(refBasesStart.charAt(i));
                ++exactMatch;
            }
            else
            {
                break;
            }
        }

        if(exactMatch == 0)
            return new HomologyData("", 0, 0, 0, overlap);

        int halfRange = exactMatch / 2;
        int exactStart = max(halfRange, 1); // round up if an odd length
        int exactEnd = exactMatch - exactStart;

        int inexactStart = halfRange;
        int inexactEnd = overlap - inexactStart;

        return new HomologyData(sb.toString(), -exactStart, exactEnd, -inexactStart, inexactEnd);
    }

    private static String getOverlapBases(
            final AlignData alignment, final Orientation orientation, final int overlap, final RefGenomeInterface refGenome)
    {
        if(orientation.isForward())
        {
            return refGenome.getBaseString(
                    alignment.refLocation().Chromosome, alignment.refLocation().end() - overlap + 1, alignment.refLocation().end());
        }
        else
        {
            return refGenome.getBaseString(
                    alignment.refLocation().Chromosome, alignment.refLocation().start(), alignment.refLocation().start() + overlap - 1);
        }
    }
}
