package com.hartwig.hmftools.sage.compare;

import java.util.Set;

import com.hartwig.hmftools.common.variant.VariantTier;

import htsjdk.variant.variantcontext.VariantContext;

public class VariantCompareData
{
    public final String Chromosome;
    public final int Position; // as per the variant definition
    public final String Ref;
    public final String Alt;

    private final VariantContext mContext;

    private final VariantTier mTier;
    private final double mQual;

    public VariantCompareData(final String chromosome, final int position, final String ref, final String alt, final VariantContext context)
    {
        Chromosome = chromosome;

        Position = position;
        Ref = ref;
        Alt = alt;

        mContext = context;
        mTier = VariantTier.fromContext(mContext);
        mQual =  mContext.getPhredScaledQual();
    }

    public double qual() { return mQual; }

    public VariantTier tier() { return mTier; }

    public Set<String> filters() { return mContext.getFilters(); }

    public String toString() { return String.format("%s:%d %s>%s tier(%s) qual(%.0f)",
            Chromosome, Position, Ref, Alt, mTier, mQual); }

    public static VariantCompareData fromContext(final VariantContext variantContext)
    {
        if(variantContext == null)
            return null;

        int variantPosition = variantContext.getStart();
        String chromosome = variantContext.getContig();

        String ref = variantContext.getReference().getBaseString();
        String alt = variantContext.getAlternateAlleles().get(0).toString();

        return new VariantCompareData(chromosome, variantPosition, ref, alt, variantContext);
    }
}
