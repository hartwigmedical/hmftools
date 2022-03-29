package com.hartwig.hmftools.common.amber;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Modifiable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public interface BaseDepthData {

    enum Base {
        G,
        A,
        T,
        C,
        N
    }

    @NotNull
    Base ref();

    @NotNull
    Base alt();

    int readDepth();

    int indelCount();

    int refSupport();

    int altSupport();

    default boolean isValid() {
        return indelCount() == 0;
    }
}
