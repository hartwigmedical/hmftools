package com.hartwig.hmftools.datamodel.purple;

import java.util.List;

import org.immutables.gson.Gson;
import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Gson.TypeAdapters
@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public interface PurpleVariant extends Variant {

    @NotNull
    PurpleCodingEffect worstCodingEffect();

    @NotNull
    PurpleTranscriptImpact canonicalImpact();

    @NotNull
    List<PurpleTranscriptImpact> otherImpacts();

    @NotNull
    Hotspot hotspot();

    boolean reported();

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

    List<PurpleVariantTranscriptImpact> variantTranscriptImpacts();

}
