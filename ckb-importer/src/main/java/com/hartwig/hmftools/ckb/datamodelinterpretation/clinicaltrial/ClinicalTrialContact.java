package com.hartwig.hmftools.ckb.datamodelinterpretation.clinicaltrial;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public abstract class ClinicalTrialContact {

    @NotNull
    public abstract String name();

    @NotNull
    public abstract String email();

    @NotNull
    public abstract String phone();

    @NotNull
    public abstract String phoneExt();

    @NotNull
    public abstract String role();
}
