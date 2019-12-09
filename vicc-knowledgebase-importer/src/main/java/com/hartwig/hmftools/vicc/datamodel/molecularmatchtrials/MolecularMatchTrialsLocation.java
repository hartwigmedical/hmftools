package com.hartwig.hmftools.vicc.datamodel.molecularmatchtrials;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public abstract class MolecularMatchTrialsLocation {

    @NotNull
    public abstract String status();

    @Nullable
    public abstract String name();

    @Nullable
    public abstract MolecularMatchTrialsContact contact();

    @Nullable
    public abstract String lastName();

    @Nullable
    public abstract String email();

    @Nullable
    public abstract String phone();

    @Nullable
    public abstract String phoneExt();

    @Nullable
    public abstract String lastNameBackup();

    @Nullable
    public abstract String emailBackup();

    @Nullable
    public abstract String phoneBackup();

    @Nullable
    public abstract String phoneExtBackup();

    @Nullable
    public abstract MolecularMatchTrialsSubLocation subLocation();

    @Nullable
    public abstract String street();

    @Nullable
    public abstract String city();

    @Nullable
    public abstract String zip();

    @Nullable
    public abstract String state();

    @Nullable
    public abstract String country();

    @Nullable
    public abstract String number();

    @Nullable
    public abstract String poBox();

    @Nullable
    public abstract String id();

    @Nullable
    public abstract String valid();

    @Nullable
    public abstract String validMessage();

    @Nullable
    public abstract String created();

    @Nullable
    public abstract String lastUpdated();

    @Nullable
    public abstract String failedGeocode();

    @Nullable
    public abstract MolecularMatchTrialsGeo geo();

}
