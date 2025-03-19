package com.hartwig.hmftools.pavereverse.gene;

import java.util.Objects;

import com.google.common.base.Preconditions;
import com.hartwig.hmftools.common.gene.ExonData;

public class AnnotatedExon
{
    public final int FirstBase;
    public final int LastBase;
    public final ExonData Exon;

    public AnnotatedExon(final int firstBase, final int lastBase, final ExonData exon)
    {
        FirstBase = firstBase;
        LastBase = lastBase;
        Exon = exon;
    }

    public boolean contains(int base)
    {
        return FirstBase <= base && LastBase >= base;
    }

    public int getAbsolutePosition(final int base)
    {
        Preconditions.checkArgument(base >= FirstBase && base <= LastBase);
        return Exon.Start + (base - FirstBase);
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
