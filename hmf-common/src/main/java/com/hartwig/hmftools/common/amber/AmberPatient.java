package com.hartwig.hmftools.common.amber;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public interface AmberPatient {

    int patientId();

    String sample();

    String otherSample();

    int matches();

    int sites();

    default double likelihood() {
        return ((double) matches()) / sites();
    }
}
