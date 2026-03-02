package com.hartwig.hmftools.amber;

import static java.lang.String.format;

import com.hartwig.hmftools.common.amber.AmberBase;
import com.hartwig.hmftools.common.amber.BaseDepthData;
import com.hartwig.hmftools.common.genome.position.GenomePosition;

public class PositionEvidence implements GenomePosition
{
    public final String Chromosome;
    public final int Position;
    public final AmberBase Ref;
    public final AmberBase Alt;

    public int ReadDepth;
    public int IndelCount;
    public int RefSupport;
    public int AltSupport;
    public int BaseQualFiltered;
    public int MapQualFiltered;

    public PositionEvidence(final String chromosome, final int position, final String ref, final String alt)
    {
        Chromosome = chromosome;
        Position = position;
        Ref = AmberBase.valueOf(ref);
        Alt = AmberBase.valueOf(alt);
        ReadDepth = 0;
        IndelCount = 0;
        RefSupport = 0;
        AltSupport = 0;
        BaseQualFiltered = 0;
        MapQualFiltered = 0;
    }

    public static PositionEvidence copy(final PositionEvidence other)
    {
        return new PositionEvidence(other.Chromosome, other.Position, other.ref(), other.alt());
    }

    public boolean isValid()
    {
        return IndelCount == 0;
    }

    public String toString()
    {
        return format("%s:%d %s>%s depth(%d) indels(%d) support(%d/%d)",
                Chromosome, Position, Ref, Alt, ReadDepth, IndelCount, RefSupport, AltSupport);
    }

    @Override
    public String chromosome()
    {
        return Chromosome;
    }

    public int position()
    {
        return Position;
    }

    public double vaf()
    {
        if(ReadDepth == 0)
        {
            return Double.NaN;
        }
        return (double) AltSupport / ReadDepth;
    }

    public String ref()
    {
        return Ref.toString();
    }

    public String alt()
    {
        return Alt.toString();
    }

    public boolean equalsRef(final char base)
    {
        return AmberBase.valueOf(String.valueOf(base)) == Ref;
    }

    public boolean equalsAlt(final char base)
    {
        return AmberBase.valueOf(String.valueOf(base)) == Alt;
    }

    public BaseDepthData toBaseDepthData()
    {
        return new BaseDepthData(
                AmberBase.valueOf(ref()),
                AmberBase.valueOf(alt()),
                ReadDepth, RefSupport, AltSupport, IndelCount);
    }
}
