package com.hartwig.hmftools.patientdb.amber;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public interface AmberSample {
    byte DO_NOT_MATCH = (byte) 0;

    String sampleId();

    byte[] entries();

}
