package com.hartwig.hmftools.patientdb.clinical.ecrf.datamodel;

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

    public static EcrfDataField of(@NotNull String patientId, @NotNull String studyOID, @NotNull String studyRepeatKey,
            @NotNull String formOID, @NotNull String formRepeatKey, @NotNull String itemGroupOID, @NotNull String itemGroupRepeatKey,
            @NotNull String itemOID, @NotNull String itemValue, @NotNull String formStatus, @NotNull String locked) {
        int studyKey = Integer.parseInt(studyRepeatKey);
        int formKey = Integer.parseInt(formRepeatKey);
        int itemGroupKey = Integer.parseInt(itemGroupRepeatKey);

        return ImmutableEcrfDataField.of(patientId,
                studyOID,
                studyKey,
                formOID,
                formKey,
                itemGroupOID,
                itemGroupKey,
                itemOID,
                itemValue,
                formStatus,
                locked);
    }
}
