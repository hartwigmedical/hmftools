package com.hartwig.hmftools.bamtools.remapper;

import java.util.Objects;
import java.util.Set;

import org.jetbrains.annotations.NotNull;

import htsjdk.samtools.SAMFlag;

public class HlaAlignmentPair implements Comparable<HlaAlignmentPair>
{
    @NotNull
    public final HlaAlignment Left;
    @NotNull
    public final HlaAlignment Right;
    private final int InterPairDistance;

    public HlaAlignmentPair(@NotNull final HlaAlignment left, @NotNull final HlaAlignment right)
    {
        this.Left = left;
        this.Right = right;
        InterPairDistance = Math.abs(Math.abs(left.Position_1Based) - Math.abs(right.Position_1Based));
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
        if (!orientationOk)
        {
            return false;
        }
        boolean sameStrand = Left.getRefId() == Right.getRefId();
        if(!sameStrand)
        {
            return false;
        }
        return InterPairDistance > 50 && InterPairDistance < 1000;
    }

    public int interPairDistance()
    {
        return InterPairDistance;
    }

    @Override
    public int compareTo(@NotNull final HlaAlignmentPair o)
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
