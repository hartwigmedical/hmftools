package com.hartwig.hmftools.lilac.variant;

import com.hartwig.hmftools.common.variant.CodingEffect;
import com.hartwig.hmftools.lilac.hla.HlaGene;

import htsjdk.variant.variantcontext.VariantContext;

public class SomaticVariant
{
    public final HlaGene Gene;
    public final String Chromosome;
    public final int Position;
    public final String Ref;
    public final String Alt;
    public final String Filter;
    public final CodingEffect CanonicalCodingEffect;

    public final VariantContext Context;

    public SomaticVariant(
            final HlaGene gene, final String chromosome, final int position, final String ref, final String alt,
            final String filter, final CodingEffect canonicalCodingEffect, final VariantContext context)
    {
        Gene = gene;
        Chromosome = chromosome;
        Position = position;
        Ref = ref;
        Alt = alt;
        Filter = filter;
        CanonicalCodingEffect = canonicalCodingEffect;
        Context = context;
    }

    @Override
    public String toString()
    {
        return String.format("gene(%s) pos(%s:%d) variant(%s>%s %s)", Gene, Chromosome, Position, Ref, Alt, CanonicalCodingEffect);
    }
}
