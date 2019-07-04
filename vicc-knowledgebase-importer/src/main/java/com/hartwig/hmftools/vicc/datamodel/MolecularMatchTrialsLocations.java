package com.hartwig.hmftools.vicc.datamodel;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public abstract class MolecularMatchTrialsLocations {

    @NotNull
    public abstract String status();

    @NotNull
    public abstract String city();

    @NotNull
    public abstract String valid();

    @NotNull
    public abstract String zip();

    @NotNull
    public abstract String created();

    @NotNull
    public abstract String country();

    @NotNull
    public abstract String id();

    @NotNull
    public abstract String lastUpdated();

    @NotNull
    public abstract MolecularMatchTrialsContact contact();

    @NotNull
    public abstract String state();

    @NotNull
    public abstract String street();

    @NotNull
    public abstract MolecularMatchTrialsLocation location();

    @NotNull
    public abstract String po_box();

    @NotNull
    public abstract String failedGeocode();

    @NotNull
    public abstract MolecularMatchTrialsGeo geo();

    @NotNull
    public abstract String validMessage();

    @NotNull
    public abstract String name();

}
