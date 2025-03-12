package com.hartwig.hmftools.pavereverse.base;

import java.util.Objects;

import com.google.common.base.Preconditions;
import com.hartwig.hmftools.common.codon.Nucleotides;
import com.hartwig.hmftools.pavereverse.util.Checks;

public class BaseSequence
{
    public final int Start;
    public final String Bases;
    public final boolean IsForwardStrand;

    public BaseSequence(int start, String bases, boolean forwardStrand)
    {
        Preconditions.checkArgument(start >= 0);
        Preconditions.checkArgument(Checks.isNucleotideSequence(bases));
        Start = start;
        Bases = bases;
        IsForwardStrand = forwardStrand;
    }

    @Override
    public String toString()
    {
        return "BaseSequence{" +
                "Start=" + Start +
                ", Bases='" + Bases + '\'' +
                ", IsForwardStrand=" + IsForwardStrand +
                '}';
    }

    @Override
    public boolean equals(final Object o)
    {
        if(o == null || getClass() != o.getClass())
        {
            return false;
        }
        final BaseSequence that = (BaseSequence) o;
        return Start == that.Start && IsForwardStrand == that.IsForwardStrand && Objects.equals(Bases, that.Bases);
    }

    @Override
    public int hashCode()
    {
        return Objects.hash(Start, Bases, IsForwardStrand);
    }

    public BaseSequence reverseComplement()
    {
        return new BaseSequence(Start, Nucleotides.reverseComplementBases(Bases), !IsForwardStrand);
    }
}
