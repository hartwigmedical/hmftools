package com.hartwig.hmftools.sage.compare;

import static com.hartwig.hmftools.common.variant.VariantVcfTags.PASS;

import java.util.Set;

import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.variant.VariantTier;

import htsjdk.variant.variantcontext.VariantContext;

public class VariantData
{
    public final String Chromosome;
    public final int Position; // as per the variant definition
    public final String Ref;
    public final String Alt;

    private final VariantContext mContext;

    private final VariantTier mTier;
    private final double mQual;

    public VariantData(final String chromosome, final int position, final String ref, final String alt, final VariantContext context)
    {
        Chromosome = chromosome;

        Position = position;
        Ref = ref;
        Alt = alt;

        mContext = context;
        mTier = VariantTier.fromContext(mContext);
        mQual =  mContext.getPhredScaledQual();
    }

    public VariantContext context() { return mContext; }
    public double qual() { return mQual; }
    public VariantTier tier() { return mTier; }
    public Set<String> filters() { return mContext.getFilters(); }

    public int allelicDepth() { return mContext.getGenotype(mContext.getGenotypes().size() - 1).getAD()[1]; }
    public int referenceDepth() { return mContext.getGenotype(mContext.getGenotypes().size() - 1).getAD()[0]; }

    public boolean isPassing()
    {
        if(mContext.getFilters().isEmpty())
            return true;

        return mContext.getFilters().stream().anyMatch(x -> x.equals(PASS));
    }

    public String toString() { return String.format("%s:%d %s>%s tier(%s) qual(%.0f)",
            Chromosome, Position, Ref, Alt, mTier, mQual); }

    public static VariantData fromContext(final VariantContext variantContext)
    {
        if(variantContext == null)
            return null;

        int variantPosition = variantContext.getStart();
        String chromosome = variantContext.getContig();

        String ref = variantContext.getReference().getBaseString();
        String alt = variantContext.getAlternateAlleles().get(0).toString();

        return new VariantData(chromosome, variantPosition, ref, alt, variantContext);
    }

    public boolean matches(final VariantData other)
    {
        return Chromosome.equals(other.Chromosome) && Position == other.Position && Ref.equals(other.Ref) && Alt.equals(other.Alt);
    }

    public static int comparePositions(final VariantData first, final VariantData second)
    {
        // return -1 if first is earlier in genomic space than second
        if(first.Chromosome.equals(second.Chromosome))
        {
            if(first.Position == second.Position)
                return 0;

            return first.Position < second.Position ? -1: 1;
        }
        else
        {
            return HumanChromosome.lowerChromosome(first.Chromosome, second.Chromosome) ? -1 : 1;
        }
    }

}
