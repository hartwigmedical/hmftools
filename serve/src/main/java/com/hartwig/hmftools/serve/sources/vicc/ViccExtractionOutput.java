package com.hartwig.hmftools.serve.sources.vicc;

import java.util.List;

import com.hartwig.hmftools.serve.actionability.fusion.ActionableFusion;
import com.hartwig.hmftools.serve.actionability.gene.ActionableGene;
import com.hartwig.hmftools.serve.actionability.hotspot.ActionableHotspot;
import com.hartwig.hmftools.serve.actionability.range.ActionableRange;
import com.hartwig.hmftools.serve.actionability.signature.ActionableSignature;
import com.hartwig.hmftools.serve.copynumber.KnownCopyNumber;
import com.hartwig.hmftools.serve.fusion.KnownFusionPair;
import com.hartwig.hmftools.serve.hotspot.KnownHotspot;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(allParameters = true,
             passAnnotations = { NotNull.class, Nullable.class })
public abstract class ViccExtractionOutput {

    @NotNull
    public abstract List<KnownHotspot> knownHotspots();

    @NotNull
    public abstract List<KnownCopyNumber> knownCopyNumbers();

    @NotNull
    public abstract List<KnownFusionPair> knownFusionPairs();

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
