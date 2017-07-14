package com.hartwig.hmftools.common.ecrf.datamodel;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(allParameters = true,
             passAnnotations = { NotNull.class, Nullable.class })
public abstract class EcrfDataField implements EcrfField {
    @NotNull
    public abstract String patientId();

    @Override
    @NotNull
    public abstract String studyEventOID();

    public abstract int studyRepeatKey();

    @Override
    @NotNull
    public abstract String formOID();

    public abstract int formRepeatKey();

    @Override
    @NotNull
    public abstract String itemGroupOID();

    public abstract int itemGroupRepeatKey();

    @Override
    @NotNull
    public abstract String itemOID();

    @NotNull
    public abstract String itemValue();

    @NotNull
    public abstract String status();

    @NotNull
    public abstract String locked();

    public static EcrfDataField of(@NotNull final String patientId, @NotNull final String studyOID, @NotNull final String studyRepeatKey,
            @NotNull final String formOID, @NotNull final String formRepeatKey, @NotNull final String itemGroupOID,
            @NotNull final String itemGroupRepeatKey, @NotNull final String itemOID, @NotNull final String itemValue,
            @NotNull final String formStatus, @NotNull final String locked) {
        final int studyKey = Integer.parseInt(studyRepeatKey);
        final int formKey = Integer.parseInt(formRepeatKey);
        final int itemGroupKey = Integer.parseInt(itemGroupRepeatKey);
        return ImmutableEcrfDataField.of(patientId, studyOID, studyKey, formOID, formKey, itemGroupOID, itemGroupKey, itemOID, itemValue,
                formStatus, locked);
    }
}
