package com.hartwig.hmftools.pave.transval;

import java.util.Objects;

import com.google.common.base.Preconditions;

import org.jetbrains.annotations.NotNull;

public class BaseSequence
{
    public final int Start;

    @NotNull
    public final String Bases;

    public BaseSequence(final int start, @NotNull final String bases)
    {
        Preconditions.checkArgument(start >= 0);
        Preconditions.checkArgument(Checks.isNucleotideSequence(bases));
        Start = start;
        this.Bases = bases;
    }

    @Override
    public boolean equals(final Object o)
    {
        if(o == null || getClass() != o.getClass())
        {
            return false;
        }
        final BaseSequence that = (BaseSequence) o;
        return Start == that.Start && Objects.equals(Bases, that.Bases);
    }

    @Override
    public int hashCode()
    {
        return Objects.hash(Start, Bases);
    }

    @Override
    public String toString()
    {
        return "BaseSequence{" +
                "Start=" + Start +
                ", bases='" + Bases + '\'' +
                '}';
    }
}
