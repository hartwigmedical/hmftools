package com.hartwig.hmftools.datamodel.cuppa;

import com.google.gson.annotations.SerializedName;
import org.immutables.gson.Gson;
import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

import java.util.List;

@Gson.TypeAdapters
@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public interface CuppaData {

    @NotNull
    List<CuppaPrediction> predictions();

    int simpleDups32To200B();

    int maxComplexSize();

    int telomericSGLs();

    @SerializedName("LINECount")
    int lineCount();
}
