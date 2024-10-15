package com.hartwig.hmftools.chord.common;

import com.google.common.annotations.VisibleForTesting;
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

    public SmallVariant(VariantContext context, String id, String chromosome, int position, String refBases, String altBases)
    {
        Context = context;

        Id = id;

        Chromosome = HumanChromosome.fromString(chromosome);
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

        Chromosome = HumanChromosome.fromString(chromosome);
        Position = position;

        RefBases = refBases;
        AltBases = altBases;
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
