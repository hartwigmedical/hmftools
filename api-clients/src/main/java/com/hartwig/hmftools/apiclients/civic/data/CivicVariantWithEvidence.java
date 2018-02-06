package com.hartwig.hmftools.apiclients.civic.data;

import static java.util.stream.Collectors.toList;

import java.util.Comparator;
import java.util.List;

import com.google.gson.annotations.SerializedName;

import org.immutables.gson.Gson;
import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Gson.TypeAdapters
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public abstract class CivicVariantWithEvidence implements CivicVariant {
    @SerializedName("evidence_items")
    public abstract List<CivicEvidenceItem> evidenceItems();

    @SerializedName("variant_aliases")
    public abstract List<String> variantAliases();

    @SerializedName("hgvs_expressions")

    public abstract List<String> hgvsExpressions();

    @Value.Lazy
    public List<CivicEvidenceItem> evidenceItemsWithDrugs() {
        return evidenceItems().stream()
                .filter(evidenceItem -> !evidenceItem.drugs().isEmpty())
                .sorted(Comparator.comparing(CivicEvidenceItem::level))
                .collect(toList());
    }

    public String summaryUrl() {
        final String URL_FORMAT = "https://civic.genome.wustl.edu/events/genes/%d/summary/variants/%d/summary#variant";
        return String.format(URL_FORMAT, geneId(), id());
    }

    @Override
    public String toString() {
        return "Civic Variant: " + name() + ": " + summaryUrl();
    }

    public boolean containsHgvsExpression(@NotNull final String hgvsExpression) {
        if (!hgvsExpression.isEmpty()) {
            for (final String variantExpression : hgvsExpressions()) {
                if (variantExpression.toLowerCase().contains(hgvsExpression.toLowerCase())) {
                    return true;
                }
            }
        }
        return false;
    }
}
