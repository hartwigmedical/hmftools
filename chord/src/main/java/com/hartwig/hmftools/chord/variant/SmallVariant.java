package com.hartwig.hmftools.chord.variant;

import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;

import htsjdk.variant.variantcontext.VariantContext;

public class SmallVariant
{
    public final VariantContext Context;

    public final String Id;

    public final HumanChromosome Chromosome;
    public final int Position;

    public final String RefBases;
    public final String AltBases;

    public SmallVariant(VariantContext context)
    {
        Context = context;

        Id = context.getID();

        Chromosome = HumanChromosome.fromString(context.getContig());
        Position = context.getStart();

        RefBases = context.getReference().getBaseString();
        AltBases = context.getAlternateAllele(0).getDisplayString();
    }

    public boolean isIndel()
    {
        return (RefBases.length()==1 && AltBases.length()>1) ||
                (AltBases.length()==1 && RefBases.length()>1);
    }

    public boolean isSnv()
    {
        return RefBases.length()==1 && AltBases.length()==1;
    }

    @Override
    public String toString()
    {
        return String.join(":", String.valueOf(Chromosome), String.valueOf(Position), RefBases, AltBases);
    }
}
