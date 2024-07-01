package com.hartwig.hmftools.sage.common;

import static java.lang.Math.abs;
import static java.lang.String.format;

import static com.hartwig.hmftools.common.utils.file.FileDelimiters.ITEM_DELIM;
import static com.hartwig.hmftools.sage.SageConstants.LONG_GERMLINE_INSERT_LENGTH;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.region.BasePosition;
import com.hartwig.hmftools.common.variant.VariantType;

public class SimpleVariant extends BasePosition
{
    public final String Ref;
    public final String Alt;
    public final VariantType Type;

    private final int mIndelLength; // number of bases inserted or deleted, negative for DELs

    public SimpleVariant(final String chromosome, final int position, final String ref, final String alt)
    {
        super(chromosome, position);
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

    public String chromosome() { return Chromosome; }
    public int position() { return Position; }

    public int positionEnd()
    {
        // the right-aligned positional end of the variant
        if(isSNV())
            return position();
        else if(isInsert())
            return position() + 1;
        else // deletes and MNVs
            return position() + Ref.length() - 1;
    }

    // convenience
    public String ref() { return Ref; }
    public String alt() { return Alt; }

    public int end() { return position() + ref().length() - 1; }

    public boolean isSNV() { return Type == VariantType.SNP; }
    public boolean isMNV() { return Type == VariantType.MNP; }
    public boolean isIndel() { return Type == VariantType.INDEL; }

    public boolean isDelete() { return ref().length() > alt().length(); }
    public boolean isInsert() { return alt().length() > ref().length(); }

    public int altLength() { return Alt.length(); }
    public int refLength() { return Ref.length(); }
    public int indelLength() { return mIndelLength; }
    public int indelLengthAbs() { return abs(mIndelLength); }

    public boolean matches(final SimpleVariant variant) { return matches(variant.Chromosome, variant.Position, variant.Ref, variant.Alt); }

    public boolean matches(final String chromosome, final int position, final String ref, final String alt)
    {
        return Chromosome.equals(chromosome) && Position == position && Ref.equals(ref) && Alt.equals(alt);
    }

    public static boolean isLongInsert(final SimpleVariant variant)
    {
        return variant.isInsert() && variant.indelLength() >= LONG_GERMLINE_INSERT_LENGTH;
    }

    public String toString() { return format("%s:%d %s>%s", Chromosome, Position, Ref, Alt); }

    private static final int VAR_ITEM_COUNT = 4;
    private static final String VAR_ITEM_DELIM = ":";

    public static List<SimpleVariant> fromConfig(final String configStr)
    {
        String[] variantStrs = configStr.split(ITEM_DELIM, -1);

        List<SimpleVariant> variants = Lists.newArrayList();

        for(String variantStr : variantStrs)
        {
            String[] items = variantStr.split(VAR_ITEM_DELIM, VAR_ITEM_COUNT);

            if(items.length != VAR_ITEM_COUNT)
                return null;

            variants.add(new SimpleVariant(items[0], Integer.parseInt(items[1]), items[2], items[3]));
        }

        return variants;
    }
}