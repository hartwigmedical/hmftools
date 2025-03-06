package com.hartwig.hmftools.pavereverse;

import java.util.Objects;

import org.jetbrains.annotations.NotNull;

public final class BaseSequenceChange
{
    @NotNull
    public final String Ref;
    @NotNull
    public final String Alt;
    @NotNull
    public final String mChromosome;
    public final int mPosition;

    public BaseSequenceChange(@NotNull final String ref, @NotNull final String alt, @NotNull final String chromosome, final int position)
    {
        this.Ref = ref;
        this.Alt = alt;
        this.mChromosome = chromosome;
        this.mPosition = position;
    }

    @NotNull
    public String toCsv()
    {
        return Ref + "," + Alt + "," + mPosition;
    }

    @Override
    public String toString()
    {
        return "BaseSequenceChange{" +
                "Ref='" + Ref + '\'' +
                ", Alt='" + Alt + '\'' +
                ", Chromosome='" + mChromosome + '\'' +
                ", Position=" + mPosition +
                '}';
    }

    @Override
    public boolean equals(final Object o)
    {
        if(o == null || getClass() != o.getClass())
        {
            return false;
        }
        final BaseSequenceChange baseSequenceChange = (BaseSequenceChange) o;
        return mPosition == baseSequenceChange.mPosition && Objects.equals(Ref, baseSequenceChange.Ref) && Objects.equals(Alt, baseSequenceChange.Alt)
                && Objects.equals(mChromosome, baseSequenceChange.mChromosome);
    }

    @Override
    public int hashCode()
    {
        return Objects.hash(Ref, Alt, mChromosome, mPosition);
    }
}
