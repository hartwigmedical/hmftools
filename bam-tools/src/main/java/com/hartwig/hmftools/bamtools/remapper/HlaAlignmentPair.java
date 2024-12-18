package com.hartwig.hmftools.bamtools.remapper;

import java.util.Objects;
import java.util.Set;

import org.broadinstitute.hellbender.utils.bwa.BwaMemAlignment;
import org.jetbrains.annotations.NotNull;

public class HlaAlignmentPair implements Comparable<HlaAlignmentPair>
{
    @NotNull
    public final HlaAlignment left;
    @NotNull
    public final HlaAlignment right;

    public HlaAlignmentPair(@NotNull final HlaAlignment left, @NotNull final HlaAlignment right)
    {
        this.left = left;
        this.right = right;
    }

    public int interPairDistance()
    {
        return Math.abs(Math.abs(left.getRefStart()) - Math.abs(right.Position));
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
        return Objects.equals(left, that.left) && Objects.equals(right, that.right);
    }

    @Override
    public int hashCode()
    {
        return Objects.hash(left, right);
    }

    @Override
    public String toString()
    {
        return "HlaAlignmentPair{" +
                "left=" + left +
                ", right=" + right +
                ", distance=" + interPairDistance() +
                '}';
    }
}
