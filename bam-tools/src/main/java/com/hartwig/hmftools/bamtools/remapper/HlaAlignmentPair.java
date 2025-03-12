package com.hartwig.hmftools.bamtools.remapper;

import static com.hartwig.hmftools.bamtools.remapper.RemapperConstants.MAX_LENGTH_FOR_PROPER_PAIR;
import static com.hartwig.hmftools.bamtools.remapper.RemapperConstants.MIN_LENGTH_FOR_PROPER_PAIR;

import java.util.Objects;
import java.util.Set;

import htsjdk.samtools.SAMFlag;

public class HlaAlignmentPair implements Comparable<HlaAlignmentPair>
{
    public final HlaAlignment Left;
    public final HlaAlignment Right;
    private final int InterPairDistance;

    public HlaAlignmentPair(final HlaAlignment left, final HlaAlignment right)
    {
        Left = left;
        Right = right;
        InterPairDistance = Math.abs(Math.abs(left.Position) - Math.abs(right.Position));
    }

    public boolean oneOfPairIsUnmapped()
    {
        return Left.isUnmapped() || Right.isUnmapped();
    }

    public boolean isConcordantPair()
    {
        final Set<SAMFlag> leftFlags = SAMFlag.getFlags(Left.getSamFlag());
        final Set<SAMFlag> rightFlags = SAMFlag.getFlags(Right.getSamFlag());
        boolean orientationOk = leftFlags.contains(SAMFlag.READ_REVERSE_STRAND) != rightFlags.contains(SAMFlag.READ_REVERSE_STRAND);
        if(!orientationOk)
        {
            return false;
        }
        boolean sameStrand = Objects.equals(Left.getRefId(), Right.getRefId());
        if(!sameStrand)
        {
            return false;
        }
        return InterPairDistance > MIN_LENGTH_FOR_PROPER_PAIR && InterPairDistance < MAX_LENGTH_FOR_PROPER_PAIR;
    }

    public int interPairDistance()
    {
        return InterPairDistance;
    }

    @Override
    public int compareTo(final HlaAlignmentPair o)
    {
        return interPairDistance() - o.interPairDistance();
    }

    @Override
    public boolean equals(final Object o)
    {
        if(o == null || getClass() != o.getClass())
        {
            return false;
        }
        final HlaAlignmentPair that = (HlaAlignmentPair) o;
        return Objects.equals(Left, that.Left) && Objects.equals(Right, that.Right);
    }

    @Override
    public int hashCode()
    {
        return Objects.hash(Left, Right);
    }

    @Override
    public String toString()
    {
        return "HlaAlignmentPair{" +
                "left=" + Left +
                ", right=" + Right +
                ", distance=" + interPairDistance() +
                '}';
    }
}
