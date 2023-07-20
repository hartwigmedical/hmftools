package com.hartwig.hmftools.orange.algo.purple;

import static com.hartwig.hmftools.common.variant.CodingEffect.MISSENSE;
import static com.hartwig.hmftools.common.variant.CodingEffect.NONE;
import static com.hartwig.hmftools.common.variant.CodingEffect.NONSENSE_OR_FRAMESHIFT;
import static com.hartwig.hmftools.common.variant.CodingEffect.SPLICE;
import static com.hartwig.hmftools.common.variant.CodingEffect.SYNONYMOUS;

import java.util.List;
import java.util.stream.Collectors;

import com.hartwig.hmftools.common.variant.CodingEffect;
import com.hartwig.hmftools.common.variant.impact.VariantEffect;

import org.jetbrains.annotations.NotNull;

public final class CodingEffectDeterminer
{
    @NotNull
    public static CodingEffect determineCodingEffect(@NotNull List<VariantEffect> variantEffects)
    {
        List<CodingEffect> simplifiedEffects = variantEffects.stream().map(CodingEffect::effect).collect(Collectors.toList());

        if(simplifiedEffects.stream().anyMatch(x -> x.equals(NONSENSE_OR_FRAMESHIFT)))
        {
            return NONSENSE_OR_FRAMESHIFT;
        }

        if(simplifiedEffects.stream().anyMatch(x -> x.equals(SPLICE)))
        {
            return SPLICE;
        }

        if(simplifiedEffects.stream().anyMatch(x -> x.equals(MISSENSE)))
        {
            return MISSENSE;
        }

        if(simplifiedEffects.stream().anyMatch(x -> x.equals(SYNONYMOUS)))
        {
            return SYNONYMOUS;
        }

        return NONE;
    }
}
