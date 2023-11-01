package com.hartwig.hmftools.common.variant;

import com.hartwig.hmftools.common.variant.impact.VariantEffect;

public enum CodingEffect
{
    NONSENSE_OR_FRAMESHIFT,
    SPLICE,
    MISSENSE,
    SYNONYMOUS,
    NONE,
    UNDEFINED;

    public static CodingEffect effect(final VariantEffect variantEffect)
    {
        switch(variantEffect)
        {
            case FRAMESHIFT:
            case STOP_GAINED:
            case STOP_LOST:
            case START_LOST:
                return NONSENSE_OR_FRAMESHIFT;

            case SPLICE_ACCEPTOR:
            case SPLICE_DONOR:
                return SPLICE;

            case MISSENSE:
            case INFRAME_DELETION:
            case INFRAME_INSERTION:
            case PHASED_INFRAME_DELETION:
            case PHASED_INFRAME_INSERTION:
            case PHASED_MISSENSE:
                return MISSENSE;

            case SYNONYMOUS:
            case PHASED_SYNONYMOUS:
                return SYNONYMOUS;
        }

        return NONE;
    }

    public static CodingEffect effect(final VariantConsequence consequence)
    {
        switch(consequence)
        {
            case FRAMESHIFT_VARIANT:
            case STOP_GAINED:
            case STOP_LOST:
            case START_LOST:
            case INITIATOR_CODON_VARIANT:
                return NONSENSE_OR_FRAMESHIFT;
            case MISSENSE_VARIANT:
            case PROTEIN_PROTEIN_CONTACT:
            case STRUCTURAL_INTERACTION_VARIANT:
            case INFRAME_DELETION:
            case INFRAME_INSERTION:
                return MISSENSE;
            case SPLICE_ACCEPTOR_VARIANT:
            case SPLICE_DONOR_VARIANT:
                return SPLICE;
            case SYNONYMOUS_VARIANT:
                return SYNONYMOUS;
        }

        return NONE;
    }

    public static boolean hasProteinImpact(final CodingEffect codingEffect)
    {
        return codingEffect == NONSENSE_OR_FRAMESHIFT || codingEffect == SPLICE || codingEffect == MISSENSE;
    }
}
