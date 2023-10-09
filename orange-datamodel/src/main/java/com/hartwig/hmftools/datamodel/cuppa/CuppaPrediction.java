package com.hartwig.hmftools.datamodel.cuppa;

import org.immutables.gson.Gson;
import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Gson.TypeAdapters
@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public interface CuppaPrediction
{
    @NotNull
    String cancerType();

    double likelihood();

    @Nullable
    Double snvPairwiseClassifier();

    @Nullable
    Double genomicPositionClassifier();

    @Nullable
    Double featureClassifier();

    @Nullable
    Double altSjCohortClassifier();

    @Nullable
    Double expressionPairwiseClassifier();
}
