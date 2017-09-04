package com.hartwig.hmftools.apiclients.civic.data;

import com.google.gson.annotations.SerializedName;
import com.hartwig.hmftools.common.variant.Variant;

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
    public abstract Long start();

    @Nullable
    public abstract Long start2();

    @Nullable
    public abstract Long stop();

    @Nullable
    public abstract Long stop2();

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

    public boolean equals(@NotNull final Variant variant) {
        final String refBases = referenceBases();
        final String altBases = variantBases();
        if (refBases == null || altBases == null) {
            return false;
        }
        return containsVariant(variant) && refBases.equals(variant.ref()) && altBases.equals(variant.alt());
    }

    private boolean checkCoordinates(@NotNull final Variant variant, final String chromosome, final Long start, final Long stop) {
        return chromosome != null && start != null && stop != null && chromosome.equals(variant.chromosome()) && start <= variant.position()
                && stop >= variant.position();

    }

    public boolean containsVariant(@NotNull final Variant variant) {
        final boolean firstCoordContainsVariant = checkCoordinates(variant, chromosome(), start(), stop());
        final boolean secondCoordContainsVariant = checkCoordinates(variant, chromosome2(), start2(), stop2());
        return firstCoordContainsVariant || secondCoordContainsVariant;
    }

    @Override
    public String toString() {
        return "chr: " + chromosome() + "[" + start() + " -> " + stop() + "]" + " / chr2: " + chromosome2() + "[" + start2() + " -> "
                + stop2() + "]";
    }
}
