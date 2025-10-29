package com.hartwig.hmftools.pavereverse;

import java.util.Objects;

public final class BaseSequenceChange
{
    public final String Ref;
    public final String Alt;
    public final String Chromosome;
    public final int Position;

    public BaseSequenceChange(String ref, String alt, String chromosome, final int position)
    {
        Ref = ref;
        Alt = alt;
        Chromosome = chromosome;
        Position = position;
    }

    public String toCsv()
    {
        return Ref + "," + Alt + "," + Position;
    }

    @Override
    public String toString()
    {
        return "BaseSequenceChange{" +
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
        final BaseSequenceChange baseSequenceChange = (BaseSequenceChange) o;
        return Position == baseSequenceChange.Position && Objects.equals(Ref, baseSequenceChange.Ref)
                && Objects.equals(Alt, baseSequenceChange.Alt)
                && Objects.equals(Chromosome, baseSequenceChange.Chromosome);
    }

    @Override
    public int hashCode()
    {
        return Objects.hash(Ref, Alt, Chromosome, Position);
    }
}
