package com.hartwig.hmftools.patientdb.clinical.ecrf.datamodel;

import com.hartwig.hmftools.patientdb.clinical.ecrf.formstatus.FormStatus;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(allParameters = true,
             passAnnotations = { NotNull.class, Nullable.class })
public abstract class ValidationFinding {

    @NotNull
    public abstract String level();

    @Nullable
    public abstract String patientIdentifier();

    @NotNull
    public abstract String message();

    @NotNull
    public abstract FormStatus formStatus();

    @NotNull
    public abstract String details();

    @NotNull
    public static ValidationFinding of(@NotNull String level, @Nullable String patientIdentifier, @NotNull String message,
            @NotNull FormStatus formStatus, @NotNull String details) {
        return ImmutableValidationFinding.of(level, patientIdentifier, message, formStatus, details);
    }

    @NotNull
    public static ValidationFinding of(@NotNull String level, @Nullable String patientIdentifier, @NotNull String message,
            @NotNull FormStatus formStatus) {
        return of(level, patientIdentifier, message, formStatus, "");
    }
}
