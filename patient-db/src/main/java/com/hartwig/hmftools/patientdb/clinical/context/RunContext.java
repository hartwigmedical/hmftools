package com.hartwig.hmftools.patientdb.clinical.context;

import org.jetbrains.annotations.NotNull;

public interface RunContext {

    @NotNull
    String runDirectory();

    @NotNull
    String setName();

    @NotNull
    String refSample();

    @NotNull
    String tumorSample();

    @NotNull
    String tumorBarcodeSample();
}
