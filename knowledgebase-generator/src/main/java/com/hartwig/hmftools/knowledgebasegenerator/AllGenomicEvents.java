package com.hartwig.hmftools.knowledgebasegenerator;

import java.util.List;
import java.util.Set;

import com.hartwig.hmftools.knowledgebasegenerator.cnv.ActionableAmplificationDeletion;
import com.hartwig.hmftools.knowledgebasegenerator.cnv.KnownAmplificationDeletion;
import com.hartwig.hmftools.knowledgebasegenerator.eventtype.EventType;
import com.hartwig.hmftools.knowledgebasegenerator.fusion.KnownFusions;
import com.hartwig.hmftools.knowledgebasegenerator.signatures.Signatures;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(allParameters = true,
             passAnnotations = { NotNull.class, Nullable.class })
public abstract class AllGenomicEvents {

    @NotNull
    public abstract List<EventType> eventType();

    @NotNull
    public abstract List<KnownAmplificationDeletion> knownAmplifications();

    @NotNull
    public abstract Set<String> uniqueAmplification();

    @NotNull
    public abstract List<KnownAmplificationDeletion> knownDeletions();

    @NotNull
    public abstract Set<String> uniqueDeletions();

    @NotNull
    public abstract List<ActionableAmplificationDeletion> actionableAmplification();

    @NotNull
    public abstract List<ActionableAmplificationDeletion> actionableDeletion();

    @NotNull
    public abstract List<Signatures> signatures();

    @NotNull
    public abstract List<KnownFusions> knownFusionPairs();

    @NotNull
    public abstract Set<String> uniqueKnownFusionPairs();

    @NotNull
    public abstract Set<String> knownFusionPromiscuousThree();

    @NotNull
    public abstract Set<String> knownFusionPromiscuousFive();

}
