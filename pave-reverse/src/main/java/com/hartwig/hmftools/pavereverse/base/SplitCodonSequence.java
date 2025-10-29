package com.hartwig.hmftools.pavereverse.base;

import static com.hartwig.hmftools.pavereverse.util.Checks.isNucleotideSequence;

import java.util.Objects;

import com.google.common.base.Preconditions;

/**
 * A sequence of bases that is subject to an insertion and deletion
 * mutation may be split across two exons. (We don't consider any
 * larger changes to be deletion insertions.)
 */
public class SplitCodonSequence
{
    private final String mLeft;
    private final String mRight;
    private final int mPositionOfChange;

    public SplitCodonSequence(String left, String right, int positionOfChange)
    {
        mPositionOfChange = positionOfChange;
        Preconditions.checkArgument(isNucleotideSequence(left));
        if(!right.isBlank())
        {
            Preconditions.checkArgument(isNucleotideSequence(right));
        }
        mLeft = left;
        mRight = right;
        Preconditions.checkArgument(completeSequence().length() % 3 == 0);
    }

    public int locationOfDeletedBases()
    {
        return mPositionOfChange;
    }

    public String retainedPrefix()
    {
        if(mRight.length() > 2)
        {
            return mLeft;
        }
        return "";
    }

    public String retainedSuffix()
    {
        if(mLeft.length() > 2)
        {
            return mRight;
        }
        return "";
    }

    boolean couldBeDeletionInsertion()
    {
        if(mRight.isBlank())
        {
            return true;
        }
        return mLeft.length() < 3 || mRight.length() < 3;
    }

    public String segmentThatIsModified()
    {
        if(mRight.isBlank())
        {
            return mLeft;
        }
        return mLeft.length() < 3 ? mRight : mLeft;
    }

    boolean spansTwoExons()
    {
        return !mRight.isBlank();
    }

    String completeSequence()
    {
        if(mRight.isBlank())
        {
            return mLeft;
        }
        return mLeft + mRight;
    }

    @Override
    public boolean equals(final Object o)
    {
        if(o == null || getClass() != o.getClass())
        {
            return false;
        }
        final SplitCodonSequence that = (SplitCodonSequence) o;
        return mPositionOfChange == that.mPositionOfChange && Objects.equals(mLeft, that.mLeft)
                && Objects.equals(mRight, that.mRight);
    }

    @Override
    public int hashCode()
    {
        return Objects.hash(mLeft, mRight, mPositionOfChange);
    }

    @Override
    public String toString()
    {
        return "SplitCodonSequence{" +
                "Left='" + mLeft + '\'' +
                ", Right='" + mRight + '\'' +
                ", PositionOfChange=" + mPositionOfChange +
                '}';
    }
}
