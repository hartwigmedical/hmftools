package com.hartwig.hmftools.common.amber;

import com.hartwig.hmftools.common.genome.position.GenomePosition;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Modifiable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public interface BaseDepth extends GenomePosition {

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
