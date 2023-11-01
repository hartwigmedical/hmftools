package com.hartwig.hmftools.datamodel.chord;

import com.google.gson.annotations.SerializedName;

import org.immutables.gson.Gson;
import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Gson.TypeAdapters
@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public interface ChordRecord
{
    @SerializedName("BRCA1Value")
    double brca1Value();

    @SerializedName("BRCA2Value")
    double brca2Value();

    double hrdValue();

    @NotNull
    ChordStatus hrStatus();

    @NotNull
    String hrdType();
}
