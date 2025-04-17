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
    CuppaMode cuppaMode();

    @NotNull
    List<CuppaPrediction> predictions();

    @NotNull
    CuppaPrediction bestPrediction();

    int simpleDups32To200B();

    int maxComplexSize();

    int telomericSGLs();

    @SerializedName("LINECount")
    int lineCount();
}