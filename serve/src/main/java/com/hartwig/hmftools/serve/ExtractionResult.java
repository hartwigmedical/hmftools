package com.hartwig.hmftools.serve;

import java.util.Set;

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
public abstract class ExtractionResult {

    @NotNull
    public abstract Set<KnownHotspot> knownHotspots();

    @NotNull
    public abstract Set<KnownCopyNumber> knownCopyNumbers();

    @NotNull
    public abstract Set<KnownFusionPair> knownFusionPairs();

    @NotNull
    public abstract Set<ActionableHotspot> actionableHotspots();

    @NotNull
    public abstract Set<ActionableRange> actionableRanges();

    @NotNull
    public abstract Set<ActionableGene> actionableGenes();

    @NotNull
    public abstract Set<ActionableFusion> actionableFusions();

    @NotNull
    public abstract Set<ActionableSignature> actionableSignatures();
}
