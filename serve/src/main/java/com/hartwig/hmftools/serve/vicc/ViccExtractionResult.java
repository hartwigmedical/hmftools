package com.hartwig.hmftools.serve.vicc;

import java.util.List;
import java.util.Map;

import com.hartwig.hmftools.common.variant.hotspot.VariantHotspot;
import com.hartwig.hmftools.serve.actionability.signature.SignatureName;
import com.hartwig.hmftools.serve.vicc.copynumber.CopyNumberAnnotation;
import com.hartwig.hmftools.serve.vicc.fusion.FusionAnnotation;
import com.hartwig.hmftools.serve.vicc.genelevel.GeneLevelAnnotation;
import com.hartwig.hmftools.serve.vicc.range.GeneRangeAnnotation;
import com.hartwig.hmftools.vicc.datamodel.Feature;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(allParameters = true,
             passAnnotations = { NotNull.class, Nullable.class })
public abstract class ViccExtractionResult {

    @NotNull
    public abstract Map<Feature, List<VariantHotspot>> hotspotsPerFeature();

    @NotNull
    public abstract Map<Feature, CopyNumberAnnotation> ampsDelsPerFeature();

    @NotNull
    public abstract Map<Feature, FusionAnnotation> fusionsPerFeature();

    @NotNull
    public abstract Map<Feature, GeneLevelAnnotation> geneLevelEventsPerFeature();

    @NotNull
    public abstract Map<Feature, List<GeneRangeAnnotation>> geneRangesPerFeature();

    @NotNull
    public abstract Map<Feature, SignatureName> signaturesPerFeature();

    @Nullable
    public abstract ActionableEvidence actionableEvidence();
}
