package com.hartwig.hmftools.datamodel.chord;

import com.google.gson.annotations.SerializedName;
import org.immutables.gson.Gson;
import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Gson.TypeAdapters
@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public abstract class ChordRecord {
    @SerializedName("BRCA1Value")
    public abstract double brca1Value();

    @SerializedName("BRCA2Value")
    public abstract double brca2Value();

    public abstract double hrdValue();

    @NotNull
    public abstract ChordStatus hrStatus();

    @NotNull
    public abstract String hrdType();
}
