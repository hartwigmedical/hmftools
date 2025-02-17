package com.hartwig.hmftools.chord.prep;

import com.google.common.annotations.VisibleForTesting;

import htsjdk.variant.variantcontext.VariantContext;

public class SmallVariant
{
    public final VariantContext Context;

    public final String Id;

    public final String Chromosome;
    public final int Position;

    public final String RefBases;
    public final String AltBases;

    public SmallVariant(VariantContext context, String id, String chromosome, int position, String refBases, String altBases)
    {
        Context = context;

        Id = id;

        Chromosome = chromosome;
        Position = position;

        RefBases = refBases;
        AltBases = altBases;
    }

    public SmallVariant(VariantContext context)
    {
        this(
                context,
                context.getID(),
                context.getContig(),
                context.getStart(),
                context.getReference().getBaseString(),
                context.getAlternateAllele(0).getDisplayString()
        );
    }

    @VisibleForTesting
    public SmallVariant(String chromosome, int position, String refBases, String altBases)
    {
        Context = null;

        Id = null;

        Chromosome = chromosome;
        Position = position;

        RefBases = refBases;
        AltBases = altBases;
    }

    @Override
    public String toString()
    {
        return String.join(":", String.valueOf(Chromosome), String.valueOf(Position), RefBases, AltBases);
    }
}
