package com.hartwig.hmftools.pave.transval;

import java.util.Objects;

import org.jetbrains.annotations.NotNull;

public final class TransvalHotspot
{
    @NotNull
    public final String Ref;
    @NotNull
    public final String Alt;
    @NotNull
    public final String Chromosome;
    public final int Position;

    public TransvalHotspot(@NotNull final String ref, @NotNull final String alt, @NotNull final String chromosome, final int position)
    {
        this.Ref = ref;
        this.Alt = alt;
        this.Chromosome = chromosome;
        this.Position = position;
    }

    @Override
    public String toString()
    {
        return "TransvalHotspot{" +
                "Ref='" + Ref + '\'' +
                ", Alt='" + Alt + '\'' +
                ", Chromosome='" + Chromosome + '\'' +
                ", Position=" + Position +
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
        return Position == transvalHotspot.Position && Objects.equals(Ref, transvalHotspot.Ref) && Objects.equals(Alt, transvalHotspot.Alt)
                && Objects.equals(Chromosome, transvalHotspot.Chromosome);
    }

    @Override
    public int hashCode()
    {
        return Objects.hash(Ref, Alt, Chromosome, Position);
    }
}
