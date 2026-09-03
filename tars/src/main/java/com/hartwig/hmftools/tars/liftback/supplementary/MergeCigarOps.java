package com.hartwig.hmftools.tars.liftback.supplementary;

import java.util.ArrayList;
import java.util.List;

import com.hartwig.hmftools.common.bam.CigarUtils;

import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;

// The CIGAR surgery a supplementary merge performs: measuring anchor runs, checking the shape either side of a
// terminal softclip, clamping a supplementary back off the primary's span, and stitching the merged spliced CIGAR.
// All of it is a pure function of a CIGAR and a few offsets, so it holds no resolver state and is shared by both
// the merge path and ref-verify.
public final class MergeCigarOps
{
    static List<CigarElement> buildMergedCigar(
            final List<CigarElement> upCigar, final List<CigarElement> downCigar,
            final int upLoss, final int downLoss, final int intronLength)
    {
        List<CigarElement> merged = new ArrayList<>(upCigar.size() + downCigar.size());
        for(int i = 0; i < upCigar.size() - 1; ++i) // upstream ops, excluding trailing S
        {
            if(i == upCigar.size() - 2 && upLoss > 0)
            {
                merged.add(new CigarElement(upCigar.get(i).getLength() - upLoss, upCigar.get(i).getOperator()));
            }
            else
            {
                merged.add(upCigar.get(i));
            }
        }
        merged.add(new CigarElement(intronLength, CigarOperator.N));
        for(int i = 1; i < downCigar.size(); ++i) // downstream ops, excluding leading S
        {
            if(i == 1 && downLoss > 0)
            {
                merged.add(new CigarElement(downCigar.get(i).getLength() - downLoss, downCigar.get(i).getOperator()));
            }
            else
            {
                merged.add(downCigar.get(i));
            }
        }
        return merged;
    }

    // Length of the M/=/X run adjacent to the terminal softclip (e.g. "5M 50S" -> 5). fromEnd scans the trailing side.
    static int matchedRun(final List<CigarElement> elements, final boolean fromEnd)
    {
        int start = fromEnd ? elements.size() - 1 : 0;
        int step = fromEnd ? -1 : 1;
        for(int i = start; i >= 0 && i < elements.size(); i += step)
        {
            CigarOperator op = elements.get(i).getOperator();
            if(op == CigarOperator.S)
                continue;
            if(op == CigarOperator.M || op == CigarOperator.EQ || op == CigarOperator.X)
            {
                return elements.get(i).getLength();
            }
            return 0;
        }
        return 0;
    }

    static boolean opAdjacentToSoftClip(final List<CigarElement> elements, final boolean leadingSide)
    {
        if(elements.size() < 2)
        {
            return false;
        }
        if(leadingSide)
        {
            if(elements.get(0).getOperator() != CigarOperator.S)
            {
                return false;
            }
            CigarOperator next = elements.get(1).getOperator();
            return next == CigarOperator.I || next == CigarOperator.D;
        }
        else
        {
            int last = elements.size() - 1;
            if(elements.get(last).getOperator() != CigarOperator.S)
            {
                return false;
            }
            CigarOperator prev = elements.get(last - 1).getOperator();
            return prev == CigarOperator.I || prev == CigarOperator.D;
        }
    }

    // ContigTranslator can expand a cross-exon M into M-N-M, leaving a supp that spans the junction itself and so
    // carries no terminal softclip to pair the primary with. Split it at an internal N and keep only the block away
    // from the primary, softclipping the rest, so it can serve as one anchor of the merge. Null when there is nothing
    // to cut: a supp starting before the primary keeps its head, one overrunning the primary's end keeps its tail.
    static ClampedSupp clampSuppToPrimaryBoundary(
            final List<CigarElement> suppCigar, final int suppStart,
            final int primaryStart, final int primaryRefEnd)
    {
        boolean keepHead = suppStart < primaryStart;
        if(!keepHead && suppStart + CigarUtils.cigarAlignedLength(suppCigar) - 1 <= primaryRefEnd)
        {
            return null;
        }

        int refCursor = suppStart;
        int readCursor = 0;
        int splitIndex = -1;
        int readAtSplit = 0;
        int refAfterSplit = -1;
        boolean lastBlockInsidePrimary = false;

        for(int i = 0; i < suppCigar.size(); ++i)
        {
            CigarOperator op = suppCigar.get(i).getOperator();
            int length = suppCigar.get(i).getLength();

            if(keepHead)
            {
                // split at the first N, but only once a later block reaches into the primary's span
                if(op == CigarOperator.N && splitIndex == -1)
                {
                    splitIndex = i;
                    readAtSplit = readCursor;
                }
                else if(splitIndex != -1 && refCursor >= primaryStart && op.isAlignment())
                {
                    return cutAt(suppCigar, suppStart, splitIndex, readAtSplit, keepHead);
                }
            }
            else
            {
                // split at the last N whose preceding block still ended inside the primary's span
                if(op.isAlignment())
                {
                    lastBlockInsidePrimary = refCursor + length - 1 <= primaryRefEnd;
                }
                else if(op == CigarOperator.N && lastBlockInsidePrimary)
                {
                    splitIndex = i;
                    readAtSplit = readCursor;
                    refAfterSplit = refCursor + length;
                }
            }

            if(op.consumesReferenceBases())
            {
                refCursor += length;
            }
            if(op.consumesReadBases())
            {
                readCursor += length;
            }
        }

        // the keep-head cut returns from inside the loop, so reaching here with one pending means no block ever
        // reached the primary
        if(keepHead || splitIndex == -1)
        {
            return null;
        }

        return cutAt(suppCigar, refAfterSplit, splitIndex, readAtSplit, keepHead);
    }

    // Drops the split N and everything on the far side of it, replaced by a softclip over the read bases that side
    // covered. The kept side keeps its own coordinates, so only a kept tail needs a new start.
    static ClampedSupp cutAt(
            final List<CigarElement> suppCigar, final int start,
            final int splitIndex, final int readAtSplit, final boolean keepHead)
    {
        List<CigarElement> trimmed = new ArrayList<>(suppCigar.size());

        if(keepHead)
        {
            trimmed.addAll(suppCigar.subList(0, splitIndex));
            int clipLength = CigarUtils.cigarBaseLength(suppCigar) - readAtSplit;
            if(clipLength > 0)
            {
                trimmed.add(new CigarElement(clipLength, CigarOperator.S));
            }
        }
        else
        {
            if(readAtSplit > 0)
            {
                trimmed.add(new CigarElement(readAtSplit, CigarOperator.S));
            }
            trimmed.addAll(suppCigar.subList(splitIndex + 1, suppCigar.size()));
        }

        return new ClampedSupp(start, trimmed);
    }

    // Moves 'shift' bases from the exon-edge M into the adjacent softclip to snap back the boundary.
    // Returns null when the edge op isn't M or trimming would leave nothing matched.
    static List<CigarElement> shiftBoundaryIntoSoftclip(
            final List<CigarElement> cigar, final int shift, final boolean rightExtend)
    {
        List<CigarElement> shifted = new ArrayList<>(cigar);
        int last = shifted.size() - 1;
        int softIdx = rightExtend ? last : 0;
        int matchIdx = rightExtend ? last - 1 : 1;
        if(matchIdx < 0 || matchIdx >= shifted.size())
        {
            return null;
        }

        CigarElement softElement = shifted.get(softIdx);
        CigarElement matchElement = shifted.get(matchIdx);
        if(matchElement.getOperator() != CigarOperator.M && matchElement.getOperator() != CigarOperator.EQ)
        {
            return null;
        }
        if(matchElement.getLength() - shift < 1)
        {
            return null;
        }

        shifted.set(matchIdx, new CigarElement(matchElement.getLength() - shift, matchElement.getOperator()));
        shifted.set(softIdx, new CigarElement(softElement.getLength() + shift, CigarOperator.S));
        return shifted;
    }

    static final class ClampedSupp
    {
        int Start;
        List<CigarElement> Cigar;

        ClampedSupp(final int start, final List<CigarElement> cigar)
        {
            Start = start;
            Cigar = cigar;
        }
    }
}
