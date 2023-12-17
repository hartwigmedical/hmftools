package com.hartwig.hmftools.datamodel.cuppa2;

import java.util.List;

import org.immutables.gson.Gson;
import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Gson.TypeAdapters
@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public interface Cuppa2Data
{
    @NotNull
    ProbabilityEntry topPrediction();

    @NotNull
    List<ProbabilityEntry> probabilities();

    @NotNull
    List<FeatureContributionEntry> featureContributions();

    @NotNull
    List<SignatureQuantileEntry> signatureQuantiles();
}
