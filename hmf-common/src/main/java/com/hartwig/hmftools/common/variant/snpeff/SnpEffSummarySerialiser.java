package com.hartwig.hmftools.common.variant.snpeff;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.variant.CodingEffect;
import com.hartwig.hmftools.common.variant.impact.VariantImpactSerialiser;

import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;

public final class SnpEffSummarySerialiser
{
    @NotNull
    public static SnpEffSummary fromDetails(@NotNull final List<String> worst, @NotNull final List<String> canonical)
    {
        final ImmutableSnpEffSummary.Builder builder = createBuilder();

        if(worst.size() == 5)
        {
            builder.worstGene(worst.get(0))
                    .worstTranscript(worst.get(1))
                    .worstEffect(VariantImpactSerialiser.readEffect(worst.get(2)))
                    .worstCodingEffect(CodingEffect.valueOf(worst.get(3)))
                    .genesAffected(Integer.parseInt(worst.get(4)));
        }

        if(canonical.size() == 6)
        {
            builder.canonicalGene(canonical.get(0))
                    .canonicalTranscript(canonical.get(1))
                    .canonicalEffect(canonical.get(2))
                    .canonicalCodingEffect(CodingEffect.valueOf(canonical.get(3)))
                    .canonicalHgvsCodingImpact(canonical.get(4))
                    .canonicalHgvsProteinImpact(canonical.get(5));
        }

        return builder.build();
    }

    @NotNull
    public static List<String> worstDetails(@NotNull final SnpEffSummary summary)
    {
        return Lists.newArrayList(summary.worstGene(),
                summary.worstTranscript(),
                VariantImpactSerialiser.writeEffect(summary.worstEffect()),
                summary.worstCodingEffect().toString(),
                String.valueOf(summary.genesAffected()));
    }

    @NotNull
    public static List<String> canonicalDetails(@NotNull final SnpEffSummary summary)
    {
        return Lists.newArrayList(summary.canonicalGene(),
                summary.canonicalTranscript(),
                VariantImpactSerialiser.writeEffect(summary.canonicalEffect()),
                summary.canonicalCodingEffect().toString(),
                summary.canonicalHgvsCodingImpact(),
                summary.canonicalHgvsProteinImpact());
    }


    @NotNull
    static ImmutableSnpEffSummary.Builder createBuilder()
    {
        return ImmutableSnpEffSummary.builder()
                .genesAffected(0)
                .worstGene(Strings.EMPTY)
                .worstEffect(Strings.EMPTY)
                .worstCodingEffect(CodingEffect.UNDEFINED)
                .worstTranscript(Strings.EMPTY)
                .canonicalGene(Strings.EMPTY)
                .canonicalEffect(Strings.EMPTY)
                .canonicalTranscript(Strings.EMPTY)
                .canonicalCodingEffect(CodingEffect.UNDEFINED)
                .canonicalHgvsCodingImpact(Strings.EMPTY)
                .canonicalHgvsProteinImpact(Strings.EMPTY);
    }
}
