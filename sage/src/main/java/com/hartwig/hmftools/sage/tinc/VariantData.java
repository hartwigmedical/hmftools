package com.hartwig.hmftools.sage.tinc;

import static java.lang.String.format;

import static com.hartwig.hmftools.common.utils.file.FileDelimiters.CSV_DELIM;
import static com.hartwig.hmftools.common.variant.SageVcfTags.AVG_BASE_QUAL;

import com.hartwig.hmftools.common.variant.GenotypeIds;
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

    public final int ReferenceABQ;
    public final int ReferenceDepth;
    public final int ReferenceAltFrags;
    public final int TumorDepth;
    public final int TumorAltFrags;

    public VariantData(
            final String chromosome, final int position, final String ref, final String alt,
            final int referenceDepth, final int referenceAltFrags, final int tumorDepth, final int tumorAltFrags)
    {
        Context = null;
        Chromosome = chromosome;
        Position = position;
        Ref = ref;
        Alt = alt;
        ReferenceABQ = 0;
        ReferenceDepth = referenceDepth;
        ReferenceAltFrags = referenceAltFrags;
        TumorDepth = tumorDepth;
        TumorAltFrags = tumorAltFrags;
    }

    public VariantData(final VariantContext variantContext, final GenotypeIds genotypeIds)
    {
        Context = variantContext;

        Position = variantContext.getStart();
        Chromosome = variantContext.getContig();

        Ref = variantContext.getReference().getBaseString();
        Alt = !variantContext.getAlternateAlleles().isEmpty() ? variantContext.getAlternateAlleles().get(0).toString() : Ref;

        Genotype refGenotype = Context.getGenotype(genotypeIds.ReferenceOrdinal);
        Genotype tumorGenotype = Context.getGenotype(genotypeIds.TumorOrdinal);
        ReferenceABQ = Integer.parseInt(refGenotype.getAnyAttribute(AVG_BASE_QUAL).toString().split(CSV_DELIM)[1]);

        ReferenceDepth = Context.getGenotype(genotypeIds.ReferenceOrdinal).getDP();
        ReferenceAltFrags = Context.getGenotype(genotypeIds.ReferenceOrdinal).getAD()[1];
        TumorDepth = Context.getGenotype(genotypeIds.TumorOrdinal).getDP();
        TumorAltFrags = Context.getGenotype(genotypeIds.TumorOrdinal).getAD()[1];
    }

    public boolean isSnv() { return Ref.length() == 1 && Alt.length() == 1; }
    public boolean isMnv() { return Ref.length() > 1 && Alt.equals(Ref); }
    public boolean isIndel() { return Alt.length() != Ref.length(); }

    public VariantTier tier() { return VariantTier.fromContext(Context); }

    public String toString()
    {
        return format("%s:%d %s>%s ref(%d/%d) tumor(%d/%d)",
                Chromosome, Position, Ref, Alt, ReferenceAltFrags, ReferenceDepth, TumorAltFrags, TumorDepth);
    }

}
