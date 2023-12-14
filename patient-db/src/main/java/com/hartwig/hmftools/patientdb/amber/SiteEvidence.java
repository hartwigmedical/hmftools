package com.hartwig.hmftools.patientdb.amber;

import static java.lang.String.format;

import static htsjdk.variant.variantcontext.Allele.NO_CALL;

import com.hartwig.hmftools.common.genome.position.GenomePosition;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;

public class SiteEvidence implements GenomePosition
{
    public final String Chromosome;
    public final int Position;
    public final String Ref;
    public final String Alt;

    public final int ReadDepth;
    public final int RefSupport;
    public final int AltSupport;

    public SiteEvidence(
            final String chromosome, final int position, final String ref, final String alt, final int readDepth,
            final int refSupport, final int altSupport)
    {
        Chromosome = chromosome;
        Position = position;
        Ref = ref;
        Alt = alt;
        ReadDepth = readDepth;
        RefSupport = refSupport;
        AltSupport = altSupport;
    }

    @Override
    public String chromosome() { return Chromosome; }
    public int position() { return Position; }

    public String toString()
    {
        return format("%s:%d %s>%s depth(%d) support(%d/%d)",
                Chromosome, Position, Ref, Alt, ReadDepth, RefSupport, AltSupport);
    }

    public static SiteEvidence fromVariantContext(final VariantContext context)
    {
        final Allele refAllele = context.getReference();
        final Genotype normal = context.getGenotype(0);

        final Allele altAllele = context.getAlternateAllele(0);

        return new SiteEvidence(
                context.getContig(), context.getStart(), refAllele.getBaseString(), altAllele.getBaseString(),
                normal.getDP(), normal.getAD()[1], normal.getAD()[0]);
    }
}
