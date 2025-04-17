package com.hartwig.hmftools.common.vcuppa;

import java.util.List;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Enclosing
@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public interface VCuppaPrediction
{
    @Value.Immutable
    interface CancerTypePrediction
    {
        String cancerType();
        double probability();
    }

    List<CancerTypePrediction> cancerTypePredictions();
}
