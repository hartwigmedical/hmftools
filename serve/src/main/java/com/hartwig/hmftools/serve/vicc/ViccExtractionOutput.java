package com.hartwig.hmftools.serve.vicc;

import java.util.List;
import java.util.Map;

import com.hartwig.hmftools.common.variant.hotspot.VariantHotspot;
import com.hartwig.hmftools.serve.actionability.fusion.ActionableFusion;
import com.hartwig.hmftools.serve.actionability.gene.ActionableGene;
import com.hartwig.hmftools.serve.actionability.hotspot.ActionableHotspot;
import com.hartwig.hmftools.serve.actionability.range.ActionableRange;
import com.hartwig.hmftools.serve.actionability.signature.ActionableSignature;
import com.hartwig.hmftools.serve.hotspot.HotspotAnnotation;
import com.hartwig.hmftools.serve.vicc.copynumber.KnownCopyNumber;
import com.hartwig.hmftools.serve.vicc.fusion.KnownFusionPair;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(allParameters = true,
             passAnnotations = { NotNull.class, Nullable.class })
public abstract class ViccExtractionOutput {

    // TODO This should be a simple list with annotated variants.
    @NotNull
    public abstract Map<VariantHotspot, HotspotAnnotation> hotspots();

    // TODO This should be a VICC-independent data class
    @NotNull
    public abstract List<KnownCopyNumber> knownAmpsDels();

    // TODO This should be a VICC-independent data class
    @NotNull
    public abstract List<KnownFusionPair> knownFusions();

    @NotNull
    public abstract List<ActionableHotspot> actionableHotspots();

    @NotNull
    public abstract List<ActionableRange> actionableRanges();

    @NotNull
    public abstract List<ActionableGene> actionableGenes();

    @NotNull
    public abstract List<ActionableFusion> actionableFusions();

    @NotNull
    public abstract List<ActionableSignature> actionableSignatures();
}
