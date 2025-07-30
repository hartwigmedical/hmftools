package com.hartwig.hmftools.pavereverse.gene;

import java.util.Objects;

import com.google.common.base.Preconditions;
import com.hartwig.hmftools.common.gene.ExonData;

public class AnnotatedExon
{
    public final int FirstBase;
    public final int LastBase;
    public final ExonData Exon;
    public final boolean IsForwardStrand;

    public AnnotatedExon(int firstBase, int lastBase, ExonData exon, boolean forwardStrand)
    {
        FirstBase = firstBase;
        LastBase = lastBase;
        Exon = exon;
        IsForwardStrand = forwardStrand;
    }

    public AnnotatedExon(final int firstBase, final int lastBase, final ExonData exon)
    {
        this(firstBase, lastBase, exon, true);
    }

    public boolean contains(int base)
    {
        return FirstBase <= base && LastBase >= base;
    }

    public int getAbsolutePositionOfBaseUpstreamOfExon(int base)
    {
        Preconditions.checkArgument(base <= 0);
        if(IsForwardStrand)
        {
            return Exon.Start + base - 1;
        }
        return Exon.End - base + 1;
    }

    public int getAbsolutePosition(final int base)
    {
        Preconditions.checkArgument(base >= FirstBase && base <= LastBase);
        final int offset = base - FirstBase;
        if(IsForwardStrand)
        {
            return Exon.Start + offset;
        }
        return Exon.End - offset;
    }

    @Override
    public boolean equals(final Object o)
    {
        if(o == null || getClass() != o.getClass())
        {
            return false;
        }
        final AnnotatedExon that = (AnnotatedExon) o;
        return FirstBase == that.FirstBase && LastBase == that.LastBase
                && Objects.equals(Exon, that.Exon);
    }

    @Override
    public int hashCode()
    {
        return Objects.hash(FirstBase, LastBase, Exon);
    }

    @Override
    public String toString()
    {
        return "AnnotatedExon{" +
                "FirstBase=" + FirstBase +
                ", LastBase=" + LastBase +
                ", Exon=" + Exon +
                '}';
    }
}
