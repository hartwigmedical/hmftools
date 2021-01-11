package com.hartwig.hmftools.ckb.drugs;

import java.util.List;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public abstract class Drugs {

    @NotNull
    public abstract String id();

    @NotNull
    public abstract String drugName();

    @NotNull
    public abstract List<String> terms();

    @NotNull
    public abstract List<String> synonyms();

    @Nullable
    public abstract String tradeName();

    @NotNull
    public abstract List<DrugDescription> drugDescriptions();

    @NotNull
    public abstract List<DrugClasses> drugClasses();

    @Nullable
    public abstract String casRegistryNum();

    @Nullable
    public abstract String nctId();

    @NotNull
    public abstract String createDate();

    @NotNull
    public abstract List<DrugsClinicalTrials> clinicalTrials();

    @NotNull
    public abstract List<DrugsEvidence> evidence();

    @NotNull
    public abstract List<DrugsTherapies> therapies();

    @NotNull
    public abstract List<DrugsGlobalApproavalStatus> globalApprovaStatus();


}
