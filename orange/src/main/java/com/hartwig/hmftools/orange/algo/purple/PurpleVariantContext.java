package com.hartwig.hmftools.orange.algo.purple;

import java.util.Collections;
import java.util.List;

import com.hartwig.hmftools.common.variant.SmallVariant;
import com.hartwig.hmftools.common.variant.impact.VariantTranscriptImpact;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public interface PurpleVariantContext extends SmallVariant
{
    @NotNull
    List<VariantTranscriptImpact> otherImpacts();

    double biallelicProbability();

    double subclonalLikelihood();

    @NotNull
    default List<String> reportableTranscriptsOrEmpty()
    {
        List<String> transcripts = reportableTranscripts();
        return transcripts != null ? transcripts : Collections.emptyList();
    }
}
