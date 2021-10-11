package com.hartwig.hmftools.patientdb.clinical.consents;

import java.util.List;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(allParameters = true,
             passAnnotations = { NotNull.class, Nullable.class })
public abstract class ConsentConfig {

    @Nullable
    public abstract String pifVersion();

    @NotNull
    public abstract List<String> cohort();

    @Nullable
    public abstract Boolean inHMF();

    @Nullable
    public abstract Boolean outsideEU();

    @Nullable
    public abstract String pif222();

    @Nullable
    public abstract List<String> pif222Values();

    @Nullable
    public abstract String pif221();

    @Nullable
    public abstract List<String> pif221Values();

    @Nullable
    public abstract String pif26HMF();

    @Nullable
    public abstract List<String> pif26HMFValues();

    @Nullable
    public abstract String pif26BUG();

    @Nullable
    public abstract List<String> pif26BUGValues();
}
