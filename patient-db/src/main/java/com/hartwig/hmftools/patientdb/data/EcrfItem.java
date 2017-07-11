package com.hartwig.hmftools.patientdb.data;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(of = "new",
             allParameters = true,
             passAnnotations = { NotNull.class, Nullable.class })
public abstract class EcrfItem {
    @NotNull
    public abstract String patientId();

    @NotNull
    public abstract String studyOID();

    public abstract int studyIdx();

    @NotNull
    public abstract String formOID();

    public abstract int formIdx();

    @NotNull
    public abstract String itemGroupOID();

    public abstract int itemGroupIdx();

    @NotNull
    public abstract String itemOID();

    @NotNull
    public abstract String itemValue();

    public abstract boolean sequenced();

    public abstract String formStatus();

    public abstract String locked();

    @NotNull
    public String sequencedString() {
        return sequenced() ? "Yes" : "No";
    }

}
