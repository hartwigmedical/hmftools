package com.hartwig.hmftools.serve.actionability.range;

import java.util.Set;

import com.hartwig.hmftools.serve.actionability.ActionableEvent;
import com.hartwig.hmftools.vicc.datamodel.Feature;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(allParameters = true,
             passAnnotations = { NotNull.class, Nullable.class })
public abstract class ActionableRange implements ActionableEvent {

    @NotNull
    public abstract String gene();

    @NotNull
    public abstract String chromosome();

    public abstract long start();

    public abstract long end();

    @NotNull
    public abstract MutationTypeFilter mutationType();

    public abstract int rangeInfo();

}
