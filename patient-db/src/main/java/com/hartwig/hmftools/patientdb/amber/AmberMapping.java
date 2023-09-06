package com.hartwig.hmftools.patientdb.amber;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public interface AmberMapping {

    String firstSample();

    String secondSample();

    int matches();

    int sites();

    default double likelihood() {
        return matches() / (double) sites();
    }

}
