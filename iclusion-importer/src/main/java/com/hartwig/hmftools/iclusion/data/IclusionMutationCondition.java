package com.hartwig.hmftools.iclusion.data;

import java.util.List;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(allParameters = true,
             passAnnotations = { NotNull.class, Nullable.class })
public abstract class IclusionMutationCondition {

    @NotNull
    public abstract List<IclusionMutation> mutations();

    @NotNull
    public abstract String logicType();
}
