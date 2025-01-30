package com.hartwig.hmftools.pave.transval;

import java.util.Objects;

import com.google.common.base.Preconditions;

import org.jetbrains.annotations.NotNull;

public final class TransvalHotspot implements Comparable<TransvalHotspot>
{
    @NotNull
    public final String Ref;
    @NotNull
    public final String Alt;
    @NotNull
    public final String mChromosome;
    public final int mPosition;

    public TransvalHotspot(@NotNull final String ref, @NotNull final String alt, @NotNull final String chromosome, final int position)
    {
        this.Ref = ref;
        this.Alt = alt;
        this.mChromosome = chromosome;
        this.mPosition = position;
    }

    @Override
    public int compareTo(@NotNull final TransvalHotspot o)
    {
        Preconditions.checkArgument(mChromosome.equals(o.mChromosome));
        int byAltLength = Alt.length() - o.Alt.length();
        if (byAltLength != 0) {
            return byAltLength;
        }
        int byRefLength = Ref.length() - o.Ref.length();
        if (byRefLength != 0) {
            return byRefLength;
        }
        int byPosition = mPosition - o.mPosition;
        if (byPosition != 0) {
            return byPosition;
        }
        return Alt.compareTo(o.Alt);
    }

    @Override
    public String toString()
    {
        return "TransvalHotspot{" +
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
        final TransvalHotspot transvalHotspot = (TransvalHotspot) o;
        return mPosition == transvalHotspot.mPosition && Objects.equals(Ref, transvalHotspot.Ref) && Objects.equals(Alt, transvalHotspot.Alt)
                && Objects.equals(mChromosome, transvalHotspot.mChromosome);
    }

    @Override
    public int hashCode()
    {
        return Objects.hash(Ref, Alt, mChromosome, mPosition);
    }
}
