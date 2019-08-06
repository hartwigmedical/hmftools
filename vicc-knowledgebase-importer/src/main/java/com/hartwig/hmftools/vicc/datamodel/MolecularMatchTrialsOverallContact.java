package com.hartwig.hmftools.vicc.datamodel;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public abstract class MolecularMatchTrialsOverallContact {

    @Nullable
    public abstract String phone();

    @Nullable
    public abstract String last_name();

    @Nullable
    public abstract String phone_ext();

    @Nullable
    public abstract String country();

    @Nullable
    public abstract String email();

    @Nullable
    public abstract String affiliation();

    @Nullable
    public abstract String city();

    @Nullable
    public abstract String name();

    @Nullable
    public abstract String zip();

    @Nullable
    public abstract String url();

    @Nullable
    public abstract String street();

    @Nullable
    public abstract String type();
}
