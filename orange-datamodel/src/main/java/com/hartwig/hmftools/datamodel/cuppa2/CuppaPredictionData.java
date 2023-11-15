package com.hartwig.hmftools.datamodel.cuppa2;

import java.util.List;

import com.google.gson.annotations.SerializedName;
import com.hartwig.hmftools.datamodel.cuppa.CuppaPrediction;

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
    List<CuppaPrediction> predictions();

    int simpleDups32To200B();

    int maxComplexSize();

    int telomericSGLs();

    @SerializedName("LINECount")
    int lineCount();
}
