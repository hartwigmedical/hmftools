package com.hartwig.hmftools.pave.compare;

import java.util.List;

import com.hartwig.hmftools.common.gene.TranscriptData;
import com.hartwig.hmftools.common.variant.VariantConsequence;
import com.hartwig.hmftools.common.variant.impact.VariantEffect;
import com.hartwig.hmftools.common.variant.snpeff.SnpEffAnnotation;
import com.hartwig.hmftools.pave.VariantTransImpact;

public final class ComparisonUtils
{
    public static SnpEffAnnotation findMatchingAnnotation(final VariantTransImpact transImpact, final List<SnpEffAnnotation> annotations)
    {
        return annotations.stream()
                .filter(x -> x.featureID() != null && x.featureID().equals(transImpact.TransData.TransName))
                .findFirst().orElse(null);
    }

    public static boolean ignoreSnpEffAnnotation(final SnpEffAnnotation annotation, final TranscriptData transData)
    {
        if(transData.CodingStart == null)
        {
            if(annotation.consequences().contains(VariantConsequence.INTRON_VARIANT))
                return true;
        }

        return false;
    }

    public static boolean effectsMatch(final VariantTransImpact transImpact, final SnpEffAnnotation annotation)
    {
        if(transImpact.effectsStr().equals(annotation.effects()))
            return true;

        if(transImpact.effects().size() != annotation.consequences().size())
            return false;

        if(transImpact.effects().stream().noneMatch(x -> annotation.consequences().stream().anyMatch(y -> effectsEqual(x, y))))
            return false;

        return true;
    }

    public static boolean effectsEqual(final VariantEffect effect, final VariantConsequence consequence)
    {
        return consequence.isParentTypeOf(effect.effect());
    }

}
