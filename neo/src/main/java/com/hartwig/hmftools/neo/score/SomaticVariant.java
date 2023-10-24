package com.hartwig.hmftools.neo.score;

import static com.hartwig.hmftools.common.neo.NeoEpitopeFile.pointMutationInfo;
import static com.hartwig.hmftools.common.variant.SageVcfTags.LIST_SEPARATOR;
import static com.hartwig.hmftools.common.variant.SageVcfTags.READ_CONTEXT_COUNT;

import com.hartwig.hmftools.common.variant.VariantReadSupport;
import com.hartwig.hmftools.common.variant.VariantType;

import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;

public class SomaticVariant
{
    public final String Chromosome;
    public final int Position;
    public final String Ref;
    public final String Alt;

    public final VariantType Type;

    public final int TumorFragments;
    public final int TumorDepth;

    public final int RnaFragments;
    public final int RnaDepth;

    public SomaticVariant(
            final String chromosome, final int position, final String ref, final String alt, final VariantType type,
            final int tumorFragments, final int tumorDepth, final int rnaFragments, final int rnaDepth)
    {
        Chromosome = chromosome;
        Position = position;
        Ref = ref;
        Alt = alt;
        Type = type;
        TumorFragments = tumorFragments;
        TumorDepth = tumorDepth;
        RnaFragments = rnaFragments;
        RnaDepth = rnaDepth;
    }

    public String variantInfo() { return pointMutationInfo(Chromosome, Position, Ref, Alt); }

    public static SomaticVariant fromContext(final VariantContext variantContext, final String sampleId, final String rnaSampleId)
    {
        int position = variantContext.getStart();
        String chromosome = variantContext.getContig();

        String ref = variantContext.getReference().getBaseString();
        String alt = !variantContext.getAlternateAlleles().isEmpty() ? variantContext.getAlternateAlleles().get(0).toString() : ref;

        Genotype tumorGenotype = variantContext.getGenotype(sampleId);
        Genotype rnaGenotype = variantContext.getGenotype(rnaSampleId);

        if(tumorGenotype == null || rnaGenotype == null)
            return null;

        VariantsCounts tumorCounts = extractVariantCounts(tumorGenotype);
        VariantsCounts rnaCounts = extractVariantCounts(rnaGenotype);

        return new SomaticVariant(
                chromosome, position, ref, alt, VariantType.type(variantContext),
                tumorCounts.VariantSupport, tumorCounts.Coverage, rnaCounts.VariantSupport, rnaCounts.Coverage);
    }

    private static VariantsCounts extractVariantCounts(final Genotype genotype)
    {
        final String[] qualCounts = genotype.getExtendedAttribute(READ_CONTEXT_COUNT, 0).toString()
                .split(LIST_SEPARATOR, -1);

        int variantSupport = 0;

        for(int i = 0; i <= VariantReadSupport.REALIGNED.ordinal(); ++i)
        {
            variantSupport += Integer.parseInt(qualCounts[i]);
        }

        return new VariantsCounts(genotype.getDP(), variantSupport);
    }

}
