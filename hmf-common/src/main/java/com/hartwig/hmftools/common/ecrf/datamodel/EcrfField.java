package com.hartwig.hmftools.common.ecrf.datamodel;

import java.util.Map;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(of = "new",
             passAnnotations = { NotNull.class, Nullable.class })
public abstract class EcrfField implements Comparable<EcrfField> {

    @Value.Parameter
    @NotNull
    public abstract String studyEventOID();

    @Value.Parameter
    @NotNull
    public abstract String formOID();

    @Value.Parameter
    @NotNull
    public abstract String itemGroupOID();

    @Value.Parameter
    @NotNull
    public abstract String itemOID();

    @Value.Parameter
    @Value.Auxiliary
    @NotNull
    public abstract String description();

    @Value.Parameter
    @Value.Auxiliary
    @NotNull
    public abstract Map<Integer, String> codeList();

    @NotNull
    public String name() {
        return EcrfFieldFunctions.name(studyEventOID(), formOID(), itemGroupOID(), itemOID());
    }

    @Override
    public int compareTo(@NotNull EcrfField other) {
        return name().compareTo(other.name());
    }
}
