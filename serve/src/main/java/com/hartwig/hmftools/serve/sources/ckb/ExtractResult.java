package com.hartwig.hmftools.serve.sources.ckb;

import java.util.List;
import java.util.Map;
import java.util.Set;

import com.hartwig.hmftools.ckb.datamodel.variant.Variant;
import com.hartwig.hmftools.common.variant.hotspot.VariantHotspot;
import com.hartwig.hmftools.serve.actionability.ActionableEvent;
import com.hartwig.hmftools.serve.extraction.codon.CodonAnnotation;
import com.hartwig.hmftools.serve.extraction.copynumber.KnownCopyNumber;
import com.hartwig.hmftools.serve.extraction.exon.ExonAnnotation;
import com.hartwig.hmftools.serve.extraction.fusion.KnownFusionPair;
import com.hartwig.hmftools.serve.extraction.gene.GeneLevelAnnotation;
import com.hartwig.hmftools.serve.extraction.signature.SignatureName;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(allParameters = true,
             passAnnotations = { NotNull.class, Nullable.class })
public abstract class ExtractResult {

    @NotNull
    public abstract Map<Variant, List<VariantHotspot>> hotspotsPerFeature();

    @NotNull
    public abstract Map<Variant, List<CodonAnnotation>> codonsPerFeature();

    @NotNull
    public abstract Map<Variant, List<ExonAnnotation>> exonsPerFeature();

    @NotNull
    public abstract Map<Variant, GeneLevelAnnotation> geneLevelEventsPerFeature();

    @NotNull
    public abstract Map<Variant, KnownCopyNumber> ampsDelsPerFeature();

    @NotNull
    public abstract Map<Variant, KnownFusionPair> fusionsPerFeature();

    @NotNull
    public abstract Map<Variant, SignatureName> signaturesPerFeature();

    @NotNull
    public abstract Set<ActionableEvent> actionableEvents();
}
