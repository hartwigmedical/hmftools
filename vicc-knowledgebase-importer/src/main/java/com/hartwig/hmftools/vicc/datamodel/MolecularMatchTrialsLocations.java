package com.hartwig.hmftools.vicc.datamodel;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public abstract class MolecularMatchTrialsLocations {

    @NotNull
    public abstract String status();

    @Nullable
    public abstract String city();

    @Nullable
    public abstract String valid();

    @Nullable
    public abstract String zip();

    @Nullable
    public abstract String created();

    @Nullable
    public abstract String country();

    @Nullable
    public abstract String id();

    @Nullable
    public abstract String lastUpdated();

    @Nullable
    public abstract MolecularMatchTrialsContact contact();

    @Nullable
    public abstract String state();

    @Nullable
    public abstract String street();

    @Nullable
    public abstract MolecularMatchTrialsLocation location();

    @Nullable
    public abstract String po_box();

    @Nullable
    public abstract String failedGeocode();

    @Nullable
    public abstract MolecularMatchTrialsGeo geo();

    @Nullable
    public abstract String validMessage();

    @Nullable
    public abstract String name();

}
