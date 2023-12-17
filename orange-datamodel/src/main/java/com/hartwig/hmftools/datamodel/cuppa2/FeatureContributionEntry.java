package com.hartwig.hmftools.datamodel.cuppa2;

import org.immutables.gson.Gson;
import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Gson.TypeAdapters
@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public interface FeatureContributionEntry
{
    String classifierName();

    String featureName();

    double featureValue();

    String cancerType();

    double featureContribution();

    int rank();

    int rankGroup();
}

