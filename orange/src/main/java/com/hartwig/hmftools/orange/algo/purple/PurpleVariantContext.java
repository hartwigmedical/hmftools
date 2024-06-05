package com.hartwig.hmftools.orange.algo.purple;

import java.util.List;

import com.hartwig.hmftools.common.genotype.GenotypeStatus;
import com.hartwig.hmftools.common.variant.AllelicDepth;
import com.hartwig.hmftools.common.variant.CodingEffect;
import com.hartwig.hmftools.common.variant.Hotspot;
import com.hartwig.hmftools.common.variant.Variant;
import com.hartwig.hmftools.common.variant.impact.VariantTranscriptImpact;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public interface PurpleVariantContext extends Variant
{
    boolean spliceRegion();

    @NotNull
    CodingEffect worstCodingEffect();

    @NotNull
    List<VariantTranscriptImpact> otherImpacts();

    @NotNull
    Hotspot hotspot();

    boolean reported();

    @Nullable
    AllelicDepth rnaDepth();

    double adjustedCopyNumber();

    double adjustedVAF();

    double minorAlleleCopyNumber();

    double variantCopyNumber();

    boolean biallelic();

    @NotNull
    GenotypeStatus genotypeStatus();

    int repeatCount();

    double subclonalLikelihood();

    @Nullable
    List<Integer> localPhaseSets();

    @NotNull
    List<String> reportableTranscripts();
}
