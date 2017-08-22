package com.hartwig.hmftools.civic.data;

import com.google.gson.annotations.SerializedName;

import org.immutables.gson.Gson;
import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Gson.TypeAdapters
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public abstract class CivicVariantMetadata {
    public abstract String name();

    public abstract int id();

    @SerializedName("evidence_items")
    protected abstract CivicEvidenceItemMetadata evidenceItems();

    public int acceptedEvidenceItems() {
        return evidenceItems().accepted();
    }

    public int rejectedEvidenceItems() {
        return evidenceItems().rejected();
    }

    public int submittedEvidenceItems() {
        return evidenceItems().submitted();
    }

    @Override
    public String toString() {
        return id() + ": " + name() + "\t[accepted: " + acceptedEvidenceItems() + ", rejected: " + rejectedEvidenceItems() + ", submitted: "
                + submittedEvidenceItems() + "].";
    }
}
