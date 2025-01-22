package com.hartwig.hmftools.pave.transval;

import static com.hartwig.hmftools.pave.transval.CodonVariant.isNucleotide;

import java.util.Objects;

import com.google.common.base.Preconditions;

import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

/**
 * A sequence of bases that is subject to an insertion and deletion
 * mutation may be split across two exons. (We don't consider any
 * larger changes to be indels.)
 * This class encapsulates the mapping of the start of a codon within
 * an exon and its possible extension into the next exon.
 */
public class SplitSequence
{
    @NotNull
    public final String left;

    @Nullable
    public final String right;

    private static boolean isNucleotideSequence(@NotNull String s)
    {
        if(s.isEmpty())
        {
            return false;
        }
        for(int i = 0; i < s.length(); i++) {
            if(!isNucleotide(s.charAt(0)))
            {
                return false;
            }
        }
        return true;
    }

    public SplitSequence(@NotNull final String left, @Nullable final String right)
    {
        Preconditions.checkArgument(isNucleotideSequence(left));
        if(right != null) {
            Preconditions.checkArgument(isNucleotideSequence(right));
        }
        this.left = left;
        this.right = right;
    }

    @Override
    public boolean equals(final Object o)
    {
        if(o == null || getClass() != o.getClass())
        {
            return false;
        }
        final SplitSequence that = (SplitSequence) o;
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
        return "SplitSequence{" +
                "left='" + left + '\'' +
                ", right='" + right + '\'' +
                '}';
    }
}
