package com.hartwig.hmftools.knowledgebasegenerator;

import java.util.List;
import java.util.Set;

import com.hartwig.hmftools.knowledgebasegenerator.cnv.ActionableAmplificationDeletion;
import com.hartwig.hmftools.knowledgebasegenerator.cnv.KnownAmplificationDeletion;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(allParameters = true,
             passAnnotations = { NotNull.class, Nullable.class })
public abstract class AllGenomicEvents {

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
}
