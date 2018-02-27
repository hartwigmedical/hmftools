package com.hartwig.hmftools.apiclients.civic.data;

import com.google.gson.annotations.SerializedName;
import com.hartwig.hmftools.common.variant.SomaticVariant;

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

    public boolean equals(@NotNull final SomaticVariant variant) {
        final String refBases = referenceBases();
        final String altBases = variantBases();
        return refBases != null && altBases != null && firstCoordinatesContainVariant(variant) && refBases.equals(variant.ref()) && altBases
                .equals(variant.alt());
    }

    private boolean coordinatesContainVariant(@NotNull final SomaticVariant variant, final String chromosome, final Long start,
            final Long stop) {
        return chromosome != null && start != null && stop != null && chromosome.equals(variant.chromosome()) && start <= variant.position()
                && stop >= variantEnd(variant);
    }

    private boolean firstCoordinatesContainVariant(@NotNull final SomaticVariant variant) {
        return coordinatesContainVariant(variant, chromosome(), start(), stop());
    }

    private boolean secondCoordinatesContainVariant(@NotNull final SomaticVariant variant) {
        return coordinatesContainVariant(variant, chromosome2(), start2(), stop2());
    }

    public boolean anyCoordinatesContainVariant(@NotNull final SomaticVariant variant) {
        return firstCoordinatesContainVariant(variant) || secondCoordinatesContainVariant(variant);
    }

    public boolean isFusion() {
        return chromosome() != null && chromosome2() != null && start() != null && start2() != null && stop() != null && stop2() != null;
    }

    public boolean isExtendedVariant(@NotNull final SomaticVariant variant) {
        final Long start = start();
        final Long stop = stop();
        return start != null && stop != null && firstCoordinatesContainVariant(variant) && start >= variant.position() - 5
                && referenceBases() == null && variantBases() == null && stop <= variantEnd(variant) + 5;
    }

    @Override
    public String toString() {
        return "chr: " + chromosome() + "[" + start() + " -> " + stop() + "]" + " / chr2: " + chromosome2() + "[" + start2() + " -> "
                + stop2() + "]";
    }

    @NotNull
    private static Long variantEnd(@NotNull final SomaticVariant variant) {
        return variant.position() + variant.ref().length() - 1;
    }

}
