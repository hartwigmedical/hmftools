package com.hartwig.hmftools.patientdb.clinical.ecrf.formstatus;

import java.util.Map;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(of = "new",
             allParameters = true,
             passAnnotations = { NotNull.class, Nullable.class })
public abstract class FormStatusModel {
    public abstract Map<FormStatusKey, FormStatus> formStatuses();

}
