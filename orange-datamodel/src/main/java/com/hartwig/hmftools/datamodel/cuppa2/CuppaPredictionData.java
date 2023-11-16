package com.hartwig.hmftools.datamodel.cuppa2;
import java.util.List;

import com.hartwig.hmftools.common.cuppa2.FeatureContributionEntry;
import com.hartwig.hmftools.common.cuppa2.ProbabilityEntry;
import com.hartwig.hmftools.common.cuppa2.SignatureQuantileEntry;

import org.immutables.gson.Gson;
import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Gson.TypeAdapters
@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public interface CuppaPredictionData
{
    @NotNull
    ProbabilityEntry topPrediction();

    @NotNull
    List<ProbabilityEntry> probs();

    @NotNull
    List<FeatureContributionEntry> featContribs();

    @NotNull
    List<SignatureQuantileEntry> sigQuantiles();
}
