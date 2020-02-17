package com.hartwig.hmftools.knowledgebasegenerator.transvar;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(allParameters = true,
             passAnnotations = { NotNull.class, Nullable.class })
public abstract class TransvarCandidateCodon {

    @NotNull
    public abstract String refCodon();

    @NotNull
    public abstract String altCodon();

}
