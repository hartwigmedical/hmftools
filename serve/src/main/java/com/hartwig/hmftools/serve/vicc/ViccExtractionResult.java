package com.hartwig.hmftools.serve.vicc;

import java.util.List;
import java.util.Map;

import com.hartwig.hmftools.common.variant.hotspot.VariantHotspot;
import com.hartwig.hmftools.serve.vicc.copynumber.KnownAmplificationDeletion;
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
    public abstract Map<Feature, KnownAmplificationDeletion> ampsDelsPerFeature();

    @NotNull
    public abstract Map<Feature, String> fusionsPerFeature();

}
