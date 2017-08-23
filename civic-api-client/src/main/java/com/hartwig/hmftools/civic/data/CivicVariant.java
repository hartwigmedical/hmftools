package com.hartwig.hmftools.civic.data;

import java.util.Comparator;
import java.util.List;
import java.util.stream.Collectors;

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

    @SerializedName("entrez_name")
    public abstract String gene();

    @SerializedName("gene_id")
    public abstract int geneId();

    @Nullable
    public abstract String name();

    @Nullable
    public abstract String description();

    @SerializedName("variant_types")
    public abstract List<CivicVariantType> variantTypes();

    public abstract CivicVariantCoordinates coordinates();

    @SerializedName("evidence_items")
    public abstract List<CivicEvidenceItem> evidenceItems();

    public List<CivicEvidenceItem> evidenceItemsWithDrugs() {
        return evidenceItems().stream()
                .filter(evidenceItem -> !evidenceItem.drugs().isEmpty())
                .sorted(Comparator.comparing(CivicEvidenceItem::level))
                .collect(Collectors.toList());
    }

    public String summaryUrl() {
        final String URL_FORMAT = "https://civic.genome.wustl.edu/events/genes/%d/summary/variants/%d/summary#variant";
        return String.format(URL_FORMAT, geneId(), id());
    }

    @Override
    public String toString() {
        return "Civic Variant: " + name() + ": " + summaryUrl();
    }
}
