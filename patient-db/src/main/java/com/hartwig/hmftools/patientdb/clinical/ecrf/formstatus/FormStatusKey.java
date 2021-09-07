package com.hartwig.hmftools.patientdb.clinical.ecrf.formstatus;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(of = "new",
             allParameters = true,
             passAnnotations = { NotNull.class, Nullable.class })
public abstract class FormStatusKey {
    @NotNull
    public abstract String patientId();

    @NotNull
    public abstract String formName();

    @NotNull
    public abstract String formSeqNum();

    @NotNull
    public abstract String studyEventName();

    @NotNull
    public abstract String studyEventSeqNum();
}
