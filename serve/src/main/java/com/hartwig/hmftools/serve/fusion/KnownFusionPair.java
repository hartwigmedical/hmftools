package com.hartwig.hmftools.serve.fusion;

import java.util.Set;

import com.hartwig.hmftools.common.serve.Knowledgebase;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(allParameters = true,
             passAnnotations = { NotNull.class, Nullable.class })
public abstract class KnownFusionPair implements FusionPair {

    @NotNull
    public abstract Set<Knowledgebase> sources();
}
