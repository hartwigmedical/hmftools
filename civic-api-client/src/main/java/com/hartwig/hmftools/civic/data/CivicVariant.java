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
public abstract class CivicVariant {

    public abstract int id();

    @Nullable
    public abstract String gene();

    @Nullable
    public abstract String name();

    @Nullable
    public abstract String description();

    @SerializedName("variant_types")
    public abstract List<CivicVariantType> variantTypes();

    public abstract CivicVariantCoordinates coordinates();

    @SerializedName("evidence_items")
    public abstract List<CivicEvidenceItem> evidenceItems();

    @Override
    public String toString() {
        return name() + "(" + id() + "): " + evidenceItems();
    }
}
