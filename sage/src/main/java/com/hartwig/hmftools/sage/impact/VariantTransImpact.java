package com.hartwig.hmftools.sage.impact;

import java.util.List;

import com.hartwig.hmftools.common.gene.TranscriptData;
import com.hartwig.hmftools.common.variant.VariantConsequence;
import com.hartwig.hmftools.common.variant.snpeff.SnpEffAnnotation;

public class VariantTransImpact
{
    public final TranscriptData TransData;
    public final VariantConsequence Consequence;
    public final String ConsequenceEffects;

    public VariantTransImpact(final TranscriptData transData, final VariantConsequence consequence)
    {
        TransData = transData;
        Consequence = consequence;
        ConsequenceEffects = consequence.description();
    }

    public VariantTransImpact(final TranscriptData transData, final String consequenceEffect)
    {
        TransData = transData;
        Consequence = VariantConsequence.fromEffect(consequenceEffect);
        ConsequenceEffects = consequenceEffect;
    }

    public String codingChange() { return ""; }
    public String proteinChange() { return ""; }

    public SnpEffAnnotation findMatchingAnnotation(final List<SnpEffAnnotation> annotations)
    {
        return annotations.stream()
                .filter(x -> x.featureID() != null && x.featureID().equals(TransData.TransName))
                .findFirst().orElse(null);
    }
}
