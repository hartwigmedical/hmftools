package com.hartwig.hmftools.pave.transval;

import static com.hartwig.hmftools.pave.transval.Checks.isNucleotideSequence;

import java.util.Objects;

import com.google.common.base.Preconditions;

import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

/**
 * A sequence of bases that is subject to an insertion and deletion
 * mutation may be split across two exons. (We don't consider any
 * larger changes to be deletion insertions.)
 * This class encapsulates the mapping of the start of a codon within
 * an exon and its possible extension into the next exon.
 */
public class SplitCodonSequence
{
    @NotNull
    private final String Left;

    @Nullable
    private final String Right;

    private final int PositionOfChange;

    public SplitCodonSequence(@NotNull final String left, @Nullable final String right, final int positionOfChange)
    {
        PositionOfChange = positionOfChange;
        Preconditions.checkArgument(isNucleotideSequence(left));
        if(right != null)
        {
            Preconditions.checkArgument(isNucleotideSequence(right));
        }
        this.Left = left;
        this.Right = right;
        Preconditions.checkArgument(completeSequence().length() % 3 == 0);
    }

    public int locationOfDeletedBases()
    {
        return PositionOfChange;
    }

    @NotNull
    public String retainedPrefix()
    {
        if(Right != null && Right.length() > 2)
        {
            return Left;
        }
        return "";
    }

    @NotNull
    public String retainedSuffix()
    {
        if (Left.length() > 2)
        {
            return Right == null ? "" : Right;
        }
        return "";
    }

    public boolean couldBeDeletionInsertion()
    {
        if(Right == null)
        {
            return true;
        }
        return Left.length() < 3 || Right.length() < 3;
    }

    @NotNull
    public String segmentThatIsModified()
    {
        if(Right == null)
        {
            return Left;
        }
        return Left.length() < 3 ? Right : Left;
    }

    public boolean spansTwoExons()
    {
        return Right != null;
    }

    @NotNull
    public String completeSequence()
    {
        if(Right == null)
        {
            return Left;
        }
        return Left + Right;
    }

    @Override
    public boolean equals(final Object o)
    {
        if(o == null || getClass() != o.getClass())
        {
            return false;
        }
        final SplitCodonSequence that = (SplitCodonSequence) o;
        return PositionOfChange == that.PositionOfChange && Objects.equals(Left, that.Left)
                && Objects.equals(Right, that.Right);
    }

    @Override
    public int hashCode()
    {
        return Objects.hash(Left, Right, PositionOfChange);
    }

    @Override
    public String toString()
    {
        return "SplitCodonicSequence{" +
                "Left='" + Left + '\'' +
                ", Right='" + Right + '\'' +
                ", PositionOfChange=" + PositionOfChange +
                '}';
    }
}
