package com.hartwig.hmftools.common.actionability.fusion;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(allParameters = true,
             passAnnotations = { NotNull.class, Nullable.class })
abstract class ActionabilityPromiscuosThree {

    @NotNull
    abstract String gene();

    @NotNull
    abstract String source();

    @NotNull
    abstract String reference();

    @NotNull
    abstract String drugsName();

    @NotNull
    abstract String drugsType();

    @NotNull
    abstract String cancerType();

    @NotNull
    abstract String level();

    @NotNull
    abstract String hmfLevel();

    @NotNull
    abstract String evidenceType();

    @NotNull
    abstract String significanceSource();

    @NotNull
    abstract String hmfResponse();
}
