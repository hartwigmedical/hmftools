package com.hartwig.hmftools.common.hla;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public interface HlaTypes
{

    @NotNull
    String status();

    @NotNull
    String typeA1();

    @NotNull
    String typeA2();

    @NotNull
    String typeB1();

    @NotNull
    String typeB2();

    @NotNull
    String typeC1();

    @NotNull
    String typeC2();

    int somaticVariants();
}
