package com.hartwig.hmftools.apiclients.civic.data;

import static java.util.stream.Collectors.groupingBy;
import static java.util.stream.Collectors.mapping;
import static java.util.stream.Collectors.toList;

import java.util.Comparator;
import java.util.List;
import java.util.Map;

import com.google.gson.annotations.SerializedName;

import org.apache.commons.lang3.tuple.ImmutablePair;
import org.apache.commons.lang3.tuple.Pair;
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

    @Value.Lazy
    public List<CivicEvidenceItem> evidenceItemsWithDrugs() {
        return evidenceItems().stream()
                .filter(evidenceItem -> !evidenceItem.drugs().isEmpty())
                .sorted(Comparator.comparing(CivicEvidenceItem::level))
                .collect(toList());
    }

    @NotNull
    @Value.Lazy
    public Map<String, Map<String, List<CivicEvidenceItem>>> groupedEvidenceItems() {
        return evidenceItemsWithDrugs().stream()
                .filter(evidenceItem -> {
                    final String evidenceDirection = evidenceItem.direction();
                    final String evidenceSignificance = evidenceItem.significance();
                    return evidenceItem.level() < 'C' && evidenceDirection != null && evidenceDirection.toLowerCase().equals("supports")
                            && evidenceSignificance != null && !evidenceSignificance.isEmpty() && !evidenceSignificance.toLowerCase()
                            .equals("n/a");
                })
                .flatMap(evidenceItem -> evidenceItem.drugs().stream().map(drug -> ImmutablePair.of(drug.name(), evidenceItem)))
                .collect(groupingBy(pair -> pair.getRight().significance(), groupingBy(Pair::getLeft, mapping(Pair::getRight, toList()))));
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
