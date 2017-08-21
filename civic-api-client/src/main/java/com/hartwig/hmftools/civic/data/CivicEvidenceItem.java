package com.hartwig.hmftools.civic.data;

import java.util.List;

import com.google.gson.annotations.SerializedName;

import org.immutables.gson.Gson;
import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Gson.TypeAdapters
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public abstract class CivicEvidenceItem {
    public abstract int id();

    public abstract String name();

    public abstract String description();

    public abstract CivicDisease disease();

    @SerializedName("evidence_level")
    public abstract Character level();

    @Nullable
    @SerializedName("clinical_significance")
    public abstract String significance();

    @Nullable
    @SerializedName("evidence_direction")
    public abstract String direction();

    //    public abstract String origin();

    public abstract List<CivicDrug> drugs();

    @Override
    public String toString() {
        return level() + ": " + direction() + " " + significance() + " to " + drugs();
    }
}
