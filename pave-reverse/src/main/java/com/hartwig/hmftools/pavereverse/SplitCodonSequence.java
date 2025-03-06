package com.hartwig.hmftools.pavereverse;

import static com.hartwig.hmftools.pavereverse.Checks.isNucleotideSequence;

import java.util.Objects;

import com.google.common.base.Preconditions;

import org.jetbrains.annotations.NotNull;

/**
 * A sequence of bases that is subject to an insertion and deletion
 * mutation may be split across two exons. (We don't consider any
 * larger changes to be deletion insertions.)
 */
class SplitCodonSequence
{
    @NotNull
    private final String Left;

    @NotNull
    private final String Right;

    private final int PositionOfChange;

    SplitCodonSequence(@NotNull final String left, @NotNull final String right, final int positionOfChange)
    {
        PositionOfChange = positionOfChange;
        Preconditions.checkArgument(isNucleotideSequence(left));
        if(!right.isBlank())
        {
            Preconditions.checkArgument(isNucleotideSequence(right));
        }
        this.Left = left;
        this.Right = right;
        Preconditions.checkArgument(completeSequence().length() % 3 == 0);
    }

    int locationOfDeletedBases()
    {
        return PositionOfChange;
    }

    @NotNull
    String retainedPrefix()
    {
        if(Right.length() > 2)
        {
            return Left;
        }
        return "";
    }

    @NotNull
    String retainedSuffix()
    {
        if (Left.length() > 2)
        {
            return Right;
        }
        return "";
    }

    boolean couldBeDeletionInsertion()
    {
        if(Right.isBlank())
        {
            return true;
        }
        return Left.length() < 3 || Right.length() < 3;
    }

    @NotNull
    String segmentThatIsModified()
    {
        if(Right.isBlank())
        {
            return Left;
        }
        return Left.length() < 3 ? Right : Left;
    }

    boolean spansTwoExons()
    {
        return !Right.isBlank();
    }

    @NotNull
    String completeSequence()
    {
        if(Right.isBlank())
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
        return "SplitCodonSequence{" +
                "Left='" + Left + '\'' +
                ", Right='" + Right + '\'' +
                ", PositionOfChange=" + PositionOfChange +
                '}';
    }
}
