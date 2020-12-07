package com.hartwig.hmftools.serve.sources.vicc;

import java.util.List;
import java.util.Map;

import com.hartwig.hmftools.common.variant.hotspot.VariantHotspot;
import com.hartwig.hmftools.serve.actionability.ActionableEvent;
import com.hartwig.hmftools.serve.actionability.signature.SignatureName;
import com.hartwig.hmftools.serve.codon.CodonAnnotation;
import com.hartwig.hmftools.serve.copynumber.KnownCopyNumber;
import com.hartwig.hmftools.serve.exon.ExonAnnotation;
import com.hartwig.hmftools.serve.fusion.KnownFusionPair;
import com.hartwig.hmftools.serve.gene.GeneLevelAnnotation;
import com.hartwig.hmftools.vicc.datamodel.Feature;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(allParameters = true,
             passAnnotations = { NotNull.class, Nullable.class })
abstract class ViccExtractionResult {

    @NotNull
    public abstract Map<Feature, List<VariantHotspot>> hotspotsPerFeature();

    @NotNull
    public abstract Map<Feature, List<CodonAnnotation>> codonsPerFeature();

    @NotNull
    public abstract Map<Feature, List<ExonAnnotation>> exonsPerFeature();

    @NotNull
    public abstract Map<Feature, GeneLevelAnnotation> geneLevelEventsPerFeature();

    @NotNull
    public abstract Map<Feature, KnownCopyNumber> ampsDelsPerFeature();

    @NotNull
    public abstract Map<Feature, KnownFusionPair> fusionsPerFeature();

    @NotNull
    public abstract Map<Feature, SignatureName> signaturesPerFeature();

    @Nullable
    public abstract ActionableEvent actionableEvent();
}
