package com.hartwig.hmftools.sage.common;

import static java.lang.String.format;

import com.hartwig.hmftools.common.genome.position.GenomePosition;
import com.hartwig.hmftools.common.variant.VariantType;

public class SimpleVariant implements GenomePosition
{
    public final String Chromosome;
    public final int Position;
    public final String Ref;
    public final String Alt;
    public final VariantType Type;

    private final int mIndelLength;

    public SimpleVariant(final String chromosome, final int position, final String ref, final String alt)
    {
        Chromosome = chromosome;
        Position = position;
        Ref = ref;
        Alt = alt;

        if(ref.length() == alt.length())
        {
            Type = ref.length() == 1 ? VariantType.SNP : VariantType.MNP;
            mIndelLength = 0;
        }
        else
        {
            Type = VariantType.INDEL;
            mIndelLength = alt().length() - ref().length();
        }
    }

    @Override
    public String chromosome() { return Chromosome; }

    @Override
    public int position() { return Position; }

    // convenience
    public String ref() { return Ref; }
    public String alt() { return Alt; }

    public int end() { return position() + ref().length() - 1; }

    public boolean isSNV() { return Type == VariantType.SNP; }
    public boolean isMNV() { return Type == VariantType.MNP; }
    public boolean isIndel() { return Type == VariantType.INDEL; }

    public boolean isDelete() { return ref().length() > alt().length(); }
    public boolean isInsert() { return alt().length() > ref().length(); }

    public int indelLength() { return mIndelLength; }

    public String toString() { return format("%s:%d %s>%s", Chromosome, Position, Ref, Alt); }

}