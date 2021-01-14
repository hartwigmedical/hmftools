package com.hartwig.hmftools.ckb.drug;

import java.util.List;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public abstract class Drug {

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
    public abstract List<DrugClass> drugClasses();

    @Nullable
    public abstract String casRegistryNum();

    @Nullable
    public abstract String nctId();

    @NotNull
    public abstract String createDate();

    @NotNull
    public abstract List<DrugClinicalTrial> clinicalTrials();

    @NotNull
    public abstract List<DrugEvidence> evidence();

    @NotNull
    public abstract List<DrugTherapy> therapies();

    @Nullable
    public abstract List<DrugGlobalApprovalStatus> globalApprovaStatus();


}
