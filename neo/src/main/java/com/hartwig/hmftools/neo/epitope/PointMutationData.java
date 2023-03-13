package com.hartwig.hmftools.neo.epitope;

import static com.hartwig.hmftools.common.variant.CodingEffect.MISSENSE;
import static com.hartwig.hmftools.common.variant.CodingEffect.NONSENSE_OR_FRAMESHIFT;

import com.hartwig.hmftools.common.variant.CodingEffect;
import com.hartwig.hmftools.common.variant.SomaticVariant;
import com.hartwig.hmftools.common.variant.VariantConsequence;

public class PointMutationData
{
    public final String Chromosome;
    public final int Position;
    public final String Ref;
    public final String Alt;
    public final String Gene;
    public final CodingEffect Effect;
    public final double VariantCopyNumber;
    public final double CopyNumber;
    public final double SubclonalLikelihood;
    public final int LocalPhaseSet;

    public PointMutationData(
            final String chromosome, final int position, final String ref, final String alt, final String gene,
            final CodingEffect effect, double variantCopyNumber, double copyNumber, double subclonalLikelihood, int localPhaseSet)
    {
        Chromosome = chromosome;
        Position = position;
        Ref = ref;
        Alt = alt;
        Gene = gene;
        Effect = effect;
        VariantCopyNumber = variantCopyNumber;
        CopyNumber = copyNumber;
        SubclonalLikelihood = subclonalLikelihood;
        LocalPhaseSet = localPhaseSet;
    }

    public static boolean isRelevantMutation(final SomaticVariant somaticVariant)
    {
        if(somaticVariant.worstCodingEffect() == NONSENSE_OR_FRAMESHIFT)
        {
            return somaticVariant.canonicalEffect().contains(VariantConsequence.FRAMESHIFT_VARIANT.parentTerm())
                    || somaticVariant.canonicalEffect().contains(VariantConsequence.STOP_LOST.parentTerm());
        }
        else if(somaticVariant.worstCodingEffect() == MISSENSE)
        {
            return true;
        }
        else
        {
            return false;
        }
    }
}
