package com.hartwig.hmftools.iclusion.datamodel;

import java.util.List;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(allParameters = true,
             passAnnotations = { NotNull.class, Nullable.class })
public abstract class IclusionTrial {

    @NotNull
    public abstract String id();

    @NotNull
    public abstract String acronym();

    @NotNull
    public abstract String title();

    @NotNull
    public abstract String eudra();

    @NotNull
    public abstract String nct();

    @NotNull
    public abstract String ipn();

    @NotNull
    public abstract String ccmo();

    @NotNull
    public abstract List<IclusionTumorLocation> tumorLocations();

    @NotNull
    public abstract List<IclusionTumorLocation> blacklistedTumorLocations();

    @NotNull
    public abstract List<IclusionMutationCondition> mutationConditions();

}
