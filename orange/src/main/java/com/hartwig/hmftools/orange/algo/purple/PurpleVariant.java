package com.hartwig.hmftools.orange.algo.purple;

import java.util.List;

import com.hartwig.hmftools.common.genotype.GenotypeStatus;
import com.hartwig.hmftools.common.purple.GermlineStatus;
import com.hartwig.hmftools.common.variant.AllelicDepth;
import com.hartwig.hmftools.common.variant.CodingEffect;
import com.hartwig.hmftools.common.variant.Hotspot;
import com.hartwig.hmftools.common.variant.VariantTier;
import com.hartwig.hmftools.common.variant.VariantType;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public abstract class PurpleVariant {

    @NotNull
    public abstract VariantType type();

    @NotNull
    public abstract String gene();

    public abstract int genesAffected();

    @NotNull
    public abstract String chromosome();

    public abstract int position();

    @NotNull
    public abstract String ref();

    @NotNull
    public abstract String alt();

    @NotNull
    public abstract CodingEffect worstCodingEffect();

    @NotNull
    public abstract PurpleTranscriptImpact canonicalImpact();

    @NotNull
    public abstract List<PurpleTranscriptImpact> otherImpacts();

    @NotNull
    public abstract Hotspot hotspot();

    public abstract boolean reported();

    public abstract boolean filtered();

    @NotNull
    public abstract String filter();

    public abstract boolean recovered();

    @NotNull
    public abstract AllelicDepth tumorDepth();

    @Nullable
    public abstract AllelicDepth rnaDepth();

    @Nullable
    public abstract AllelicDepth referenceDepth();

    public abstract double adjustedCopyNumber();

    public abstract double adjustedVAF();

    public abstract double minorAlleleCopyNumber();

    public abstract double variantCopyNumber();

    public abstract boolean biallelic();

    @NotNull
    public abstract GenotypeStatus genotypeStatus();

    @NotNull
    public abstract GermlineStatus germlineStatus();

    @NotNull
    public abstract String trinucleotideContext();

    public abstract double mappability();

    @NotNull
    public abstract String microhomology();

    @NotNull
    public abstract String repeatSequence();

    public abstract int repeatCount();

    @NotNull
    public abstract String kataegis();

    @NotNull
    public abstract VariantTier tier();

    public abstract double subclonalLikelihood();

    @Nullable
    public abstract List<Integer> localPhaseSets();
}
