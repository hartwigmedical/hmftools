package com.hartwig.hmftools.pave.transval;

import java.util.Objects;

import com.google.common.base.Preconditions;

import org.jetbrains.annotations.NotNull;

public class CodonVariant implements Comparable<CodonVariant>
{
    @NotNull
    public final String referenceCodon;
    @NotNull
    public final String alternateCodon;

    private static boolean isCodon(@NotNull String s)
    {
        if(s.length() != 3)
        {
            return false;
        }
        return isNucleotide(s.charAt(0)) && isNucleotide(s.charAt(1)) && isNucleotide(s.charAt(2));
    }

    private static boolean isNucleotide(char c)
    {
        return c == 'A' || c == 'C' || c == 'G' || c == 'T';
    }

    public CodonVariant(@NotNull final String referenceCodon, @NotNull final String alternateCodon)
    {
        Preconditions.checkArgument(isCodon(referenceCodon));
        Preconditions.checkArgument(isCodon(alternateCodon));
        this.referenceCodon = referenceCodon;
        this.alternateCodon = alternateCodon;
    }

    @Override
    public int compareTo(@NotNull final CodonVariant o)
    {
        if(!referenceCodon.equals(o.referenceCodon))
        {
            throw new IllegalArgumentException(String.format("%s != %s", referenceCodon, o.referenceCodon));
        }
        int byDistance = editDistance() - o.editDistance();
        if(byDistance != 0)
        {
            return byDistance;
        }
        return alternateCodon.compareTo(o.alternateCodon);
    }

    public int editDistance()
    {
        int result = 0;
        for(int i = 0; i < 3; ++i)
        {
            if(referenceCodon.charAt(i) != alternateCodon.charAt(i))
            {
                result += 1;
            }
        }
        return result;
    }

    @Override
    public boolean equals(final Object o)
    {
        if(o == null || getClass() != o.getClass())
        {
            return false;
        }
        final CodonVariant that = (CodonVariant) o;
        return Objects.equals(referenceCodon, that.referenceCodon) && Objects.equals(alternateCodon, that.alternateCodon);
    }

    @Override
    public int hashCode()
    {
        return Objects.hash(referenceCodon, alternateCodon);
    }

    @Override
    public String toString()
    {
        return "CodonVariant{" +
                "referenceCodon='" + referenceCodon + '\'' +
                ", alternateCodon='" + alternateCodon + '\'' +
                '}';
    }
}
