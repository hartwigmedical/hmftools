package com.hartwig.hmftools.sage.tinc;

import static java.lang.String.format;

import com.hartwig.hmftools.common.variant.VariantTier;
import com.hartwig.hmftools.common.variant.VariantType;

import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;

public class VariantData
{
    public final VariantContext Context;
    public final String Chromosome;
    public final int Position;
    public final String Ref;
    public final String Alt;

    public VariantData(final VariantContext variantContext)
    {
        Context = variantContext;

        Position = variantContext.getStart();
        Chromosome = variantContext.getContig();

        Ref = variantContext.getReference().getBaseString();
        Alt = !variantContext.getAlternateAlleles().isEmpty() ? variantContext.getAlternateAlleles().get(0).toString() : Ref;
    }

    public boolean isSnv() { return Ref.length() == 1 && Alt.length() == 1; }
    public boolean isMnv() { return Ref.length() > 1 && Alt.equals(Ref); }
    public boolean isIndel() { return Alt.length() != Ref.length(); }

    public VariantTier tier() { return VariantTier.fromContext(Context); }

    /*
    public int getGenotypeIntValue(final Genotype genotype, final String vcfTag)
    {

    }
    */

    public String toString()
    {
        return format("%s:%d %s>%s", Chromosome, Position, Ref, Alt);
    }

}
