package com.hartwig.hmftools.datamodel.purple;

import java.util.List;

import org.immutables.gson.Gson;
import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Gson.TypeAdapters
@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public interface PurpleVariant
{
    @NotNull
    PurpleVariantType type();

    @NotNull
    String gene();

    @NotNull
    String chromosome();

    int position();

    @NotNull
    String ref();

    @NotNull
    String alt();

    @NotNull
    PurpleCodingEffect worstCodingEffect();

    @NotNull
    PurpleTranscriptImpact canonicalImpact();

    @NotNull
    List<PurpleTranscriptImpact> otherImpacts();

    @NotNull
    HotspotType hotspot();

    @NotNull
    PurpleAllelicDepth tumorDepth();

    @Nullable
    PurpleAllelicDepth rnaDepth();

    double adjustedCopyNumber();

    double adjustedVAF();

    double minorAlleleCopyNumber();

    double variantCopyNumber();

    boolean biallelic();

    @NotNull
    PurpleGenotypeStatus genotypeStatus();

    int repeatCount();

    double subclonalLikelihood();

    @Nullable
    List<Integer> localPhaseSets();

    @Value.Derived
    default boolean reported()
    {
        if(canonicalImpact().reported())
        {
            return true;
        }
        else
        {
            return otherImpacts().stream().anyMatch(PurpleTranscriptImpact::reported);
        }
    }
}
