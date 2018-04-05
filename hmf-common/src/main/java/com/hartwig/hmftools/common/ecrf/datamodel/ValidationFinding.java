package com.hartwig.hmftools.common.ecrf.datamodel;

import com.hartwig.hmftools.common.ecrf.formstatus.FormStatus;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(allParameters = true,
             passAnnotations = { NotNull.class, Nullable.class })
public abstract class ValidationFinding {
    @NotNull
    public abstract String level();

    @NotNull
    public abstract String patientId();

    @NotNull
    public abstract String ecrfItem();

    @NotNull
    public abstract String message();

    @NotNull
    public abstract FormStatus formStatus();

    @NotNull
    public abstract String details();

    @NotNull
    public static ValidationFinding of(@NotNull String level, @NotNull String patientId, @NotNull String ecrfItem, @NotNull String message,
            @NotNull FormStatus formStatus, @NotNull String details) {
        return ImmutableValidationFinding.of(level, patientId, ecrfItem, message, formStatus, details);
    }

    @NotNull
    public static ValidationFinding of(@NotNull String level, @NotNull String patientId, @NotNull String ecrfItem, @NotNull String message,
            @NotNull FormStatus formStatus) {
        return of(level, patientId, ecrfItem, message, formStatus, "");
    }
}
