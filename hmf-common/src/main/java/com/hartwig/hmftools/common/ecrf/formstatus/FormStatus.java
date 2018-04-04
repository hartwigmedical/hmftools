package com.hartwig.hmftools.common.ecrf.formstatus;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(of = "new",
             allParameters = true,
             passAnnotations = { NotNull.class, Nullable.class })
public abstract class FormStatus {

    @NotNull
    public abstract FormStatusState state();

    public abstract boolean locked();

    @NotNull
    public static FormStatus unknown() {
        return ImmutableFormStatus.builder().state(FormStatusState.UNKNOWN).locked(false).build();
    }
}
