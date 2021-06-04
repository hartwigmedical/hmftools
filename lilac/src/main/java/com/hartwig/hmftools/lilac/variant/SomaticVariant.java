package com.hartwig.hmftools.lilac.variant;

import com.hartwig.hmftools.common.variant.CodingEffect;

import org.apache.commons.math3.distribution.PoissonDistribution;

import htsjdk.variant.variantcontext.VariantContext;

public class SomaticVariant
{
    public final String Gene;
    public final String Chromosome;
    public final int Position;
    public final String Ref;
    public final String Alt;
    public final String Filter;
    public final CodingEffect CanonicalCodingEffect;

    public final VariantContext Context;

    public SomaticVariant(
            final String gene, final String chromosome, final int position, final String ref, final String alt,
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

    public String toString()
    {
        return String.format("gene(%s) pos(%s:%d) variant(%s>%s %s)", Gene, Chromosome, Position, Ref, Alt, CanonicalCodingEffect);
    }
}
