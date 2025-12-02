package com.hartwig.hmftools.datamodel.finding;

import java.util.ArrayList;
import java.util.List;
import java.util.Set;

import com.hartwig.hmftools.datamodel.purple.PurpleTranscriptImpact;
import com.hartwig.hmftools.datamodel.purple.PurpleVariant;
import com.hartwig.hmftools.datamodel.purple.PurpleVariantEffect;

import org.jetbrains.annotations.NotNull;

public final class VariantDedup
{
    private static final Set<PurpleVariantEffect> PHASED_EFFECTS =
            Set.of(PurpleVariantEffect.PHASED_INFRAME_DELETION, PurpleVariantEffect.PHASED_INFRAME_INSERTION,
                    PurpleVariantEffect.PHASED_MISSENSE, PurpleVariantEffect.PHASED_SYNONYMOUS);

    @NotNull
    public static List<PurpleVariant> apply(@NotNull List<PurpleVariant> variants)
    {
        List<PurpleVariant> filtered = new ArrayList<>();
        for(PurpleVariant variant : variants)
        {
            if(!hasCanonicalPhasedEffect(variant) || !hasSameEffectWithHigherVCN(variants, variant))
            {
                filtered.add(variant);
            }
        }
        return filtered;
    }

    private static boolean hasCanonicalPhasedEffect(@NotNull PurpleVariant variant)
    {
        for(PurpleVariantEffect effect : variant.canonicalImpact().effects())
        {
            if(PHASED_EFFECTS.contains(effect))
            {
                return true;
            }
        }
        return false;
    }

    private static boolean hasSameEffectWithHigherVCN(@NotNull List<PurpleVariant> variants, @NotNull PurpleVariant variantToMatch)
    {
        // We assume that variants with same effect have unique hgvs coding impact.
        Double minVariantCopyNumber = null;
        String uniqueHgvsCodingImpact = null;
        PurpleTranscriptImpact variantImpactToMatch = variantToMatch.canonicalImpact();

        for(PurpleVariant variant : variants)
        {
            PurpleTranscriptImpact variantImpact = variant.canonicalImpact();
            if(variantImpact.effects().equals(variantImpactToMatch.effects()) && variant.gene().equals(variantToMatch.gene())
                    && variantImpact.hgvsProteinImpact().equals(variantImpactToMatch.hgvsProteinImpact()))
            {
                if(minVariantCopyNumber == null || Doubles.lessThan(variant.variantCopyNumber(), minVariantCopyNumber))
                {
                    minVariantCopyNumber = variant.variantCopyNumber();
                    uniqueHgvsCodingImpact = variantImpact.hgvsCodingImpact();
                }
                else if(Doubles.equal(variant.variantCopyNumber(), minVariantCopyNumber))
                {
                    uniqueHgvsCodingImpact = variantImpact.hgvsCodingImpact().compareTo(uniqueHgvsCodingImpact) > 0
                            ? variantImpact.hgvsCodingImpact()
                            : uniqueHgvsCodingImpact;
                }
            }
        }

        boolean matchesMinVariantCopyNumber = Doubles.equal(variantToMatch.variantCopyNumber(), minVariantCopyNumber);
        boolean matchesBestHgvsCodingImpact = variantImpactToMatch.hgvsCodingImpact().equals(uniqueHgvsCodingImpact);
        return !(matchesMinVariantCopyNumber && matchesBestHgvsCodingImpact);
    }
}
