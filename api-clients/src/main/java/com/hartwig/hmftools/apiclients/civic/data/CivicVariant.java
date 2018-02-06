package com.hartwig.hmftools.apiclients.civic.data;

import java.util.List;

import com.google.gson.annotations.SerializedName;

import org.immutables.gson.Gson;
import org.jetbrains.annotations.Nullable;

@Gson.TypeAdapters
public interface CivicVariant {
    int id();

    @SerializedName("entrez_name")
    String gene();

    @SerializedName("entrez_id")
    String entrezId();

    @Nullable
    String name();

    @Nullable
    String description();

    @SerializedName("gene_id")
    int geneId();

    String type();

    @SerializedName("variant_types")
    List<CivicVariantType> variantTypes();

    CivicVariantCoordinates coordinates();
}
