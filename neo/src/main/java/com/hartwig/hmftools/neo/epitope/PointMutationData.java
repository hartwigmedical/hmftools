package com.hartwig.hmftools.neo.epitope;

import static com.hartwig.hmftools.common.variant.CodingEffect.MISSENSE;
import static com.hartwig.hmftools.common.variant.CodingEffect.NONE;
import static com.hartwig.hmftools.common.variant.CodingEffect.NONSENSE_OR_FRAMESHIFT;

import com.hartwig.hmftools.common.variant.CodingEffect;
import com.hartwig.hmftools.common.variant.VariantConsequence;
import com.hartwig.hmftools.common.variant.impact.VariantImpact;
import com.hartwig.hmftools.common.variant.impact.VariantTranscriptImpact;

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

    public PointMutationData(
            final String chromosome, final int position, final String ref, final String alt, final String gene,
            final CodingEffect effect, double variantCopyNumber, double copyNumber, double subclonalLikelihood)
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
    }

    public static boolean isRelevantMutation(final VariantImpact impact)
    {
        if(impact.WorstCodingEffect == NONSENSE_OR_FRAMESHIFT)
        {
            return impact.CanonicalEffect.contains(VariantConsequence.FRAMESHIFT_VARIANT.parentTerm())
                    || impact.CanonicalEffect.contains(VariantConsequence.STOP_LOST.parentTerm());
        }
        else if(impact.WorstCodingEffect == MISSENSE)
        {
            return true;
        }
        else
        {
            return false;
        }
    }

    public static CodingEffect checkVariantEffects(final VariantTranscriptImpact transcriptImpact)
    {
        if(transcriptImpact.Effects.contains(VariantConsequence.FRAMESHIFT_VARIANT.parentTerm()))
            return NONSENSE_OR_FRAMESHIFT;

        if(transcriptImpact.Effects.contains(VariantConsequence.STOP_LOST.parentTerm()))
            return NONSENSE_OR_FRAMESHIFT;

        if(transcriptImpact.Effects.contains(VariantConsequence.MISSENSE_VARIANT.parentTerm()))
            return MISSENSE;

        return NONE;
    }
}
