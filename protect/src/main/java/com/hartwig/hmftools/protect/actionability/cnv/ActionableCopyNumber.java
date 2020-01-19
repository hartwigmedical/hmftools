package com.hartwig.hmftools.protect.actionability.cnv;

import com.hartwig.hmftools.protect.actionability.Actionable;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(allParameters = true,
             passAnnotations = { NotNull.class, Nullable.class })
public abstract class ActionableCopyNumber implements Actionable {

    @NotNull
    public abstract String gene();

    @NotNull
    public abstract CopyNumberType type();

}
