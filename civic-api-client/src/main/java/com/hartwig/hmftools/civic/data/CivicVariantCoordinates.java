package com.hartwig.hmftools.civic.data;

import com.google.gson.annotations.SerializedName;

import org.immutables.gson.Gson;
import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Gson.TypeAdapters
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public abstract class CivicVariantCoordinates {
    @Nullable
    public abstract String chromosome();

    @Nullable
    public abstract String chromosome2();

    @Nullable
    public abstract Integer start();

    @Nullable
    public abstract Integer start2();

    @Nullable
    public abstract Integer stop();

    @Nullable
    public abstract Integer stop2();

    @Nullable
    @SerializedName("representative_transcript")
    public abstract String representativeTranscript();

    @Nullable
    @SerializedName("representative_transcript2")
    public abstract String representativeTranscript2();

    @Nullable
    @SerializedName("reference_bases")
    public abstract String referenceBases();

    @Nullable
    @SerializedName("variant_bases")
    public abstract String variantBases();

    @Nullable
    @SerializedName("ensembl_version")
    public abstract String ensemblVersion();

    @Nullable
    @SerializedName("reference_build")
    public abstract String referenceBuild();
}
