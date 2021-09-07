package com.hartwig.hmftools.patientdb.clinical.ecrf.datamodel;

import java.util.Map;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(allParameters = true,
             of = "new",
             passAnnotations = { NotNull.class, Nullable.class })
public abstract class EcrfDatamodelField implements Comparable<EcrfDatamodelField>, EcrfField {

    @Override
    @NotNull
    public abstract String studyEventOID();

    @Override
    @NotNull
    public abstract String formOID();

    @Override
    @NotNull
    public abstract String itemGroupOID();

    @Override
    @NotNull
    public abstract String itemOID();

    @Value.Auxiliary
    @NotNull
    public abstract String description();

    @Value.Auxiliary
    @NotNull
    public abstract Map<Integer, String> codeList();

    @Override
    public int compareTo(@NotNull EcrfDatamodelField other) {
        return name().compareTo(other.name());
    }
}
