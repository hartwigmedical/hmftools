package com.hartwig.hmftools.datamodel.purple;

import org.immutables.gson.Gson;
import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

import java.util.List;

@Gson.TypeAdapters
@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public abstract class PurpleVariant implements Variant {

    @NotNull
    public abstract PurpleVariantType type();

    @NotNull
    public abstract String gene();

    @NotNull
    public abstract String chromosome();

    public abstract int position();

    @NotNull
    public abstract String ref();

    @NotNull
    public abstract String alt();

    @NotNull
    public abstract PurpleCodingEffect worstCodingEffect();

    @NotNull
    public abstract PurpleTranscriptImpact canonicalImpact();

    @NotNull
    public abstract List<PurpleTranscriptImpact> otherImpacts();

    @NotNull
    public abstract Hotspot hotspot();

    public abstract boolean reported();

    @NotNull
    public abstract PurpleAllelicDepth tumorDepth();

    @Nullable
    public abstract PurpleAllelicDepth rnaDepth();

    public abstract double adjustedCopyNumber();

    public abstract double adjustedVAF();

    public abstract double minorAlleleCopyNumber();

    public abstract double variantCopyNumber();

    public abstract boolean biallelic();

    @NotNull
    public abstract PurpleGenotypeStatus genotypeStatus();

    public abstract int repeatCount();

    public abstract double subclonalLikelihood();

    @Nullable
    public abstract List<Integer> localPhaseSets();
}
