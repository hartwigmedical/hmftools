package com.hartwig.hmftools.datamodel.cuppa;

import java.util.List;

import com.google.gson.annotations.SerializedName;

import org.immutables.gson.Gson;
import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Gson.TypeAdapters
@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public interface CuppaData
{
    @NotNull
    List<CuppaPrediction> predictions();

    @Nullable
    CuppaPrediction bestPrediction();

    @Nullable
    Integer simpleDups32To200B();

    @Nullable
    Integer maxComplexSize();

    @Nullable
    Integer telomericSGLs();

    @Nullable
    @SerializedName("LINECount")
    Integer lineCount();
}
