package com.hartwig.hmftools.serve.sources.iclusion;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.iclusion.datamodel.IclusionMutation;
import com.hartwig.hmftools.iclusion.datamodel.IclusionMutationCondition;
import com.hartwig.hmftools.iclusion.datamodel.IclusionMutationLogicType;
import com.hartwig.hmftools.iclusion.datamodel.IclusionTrial;
import com.hartwig.hmftools.iclusion.datamodel.IclusionTumorLocation;
import com.hartwig.hmftools.iclusion.datamodel.ImmutableIclusionMutation;
import com.hartwig.hmftools.iclusion.datamodel.ImmutableIclusionMutationCondition;
import com.hartwig.hmftools.iclusion.datamodel.ImmutableIclusionTrial;

import org.jetbrains.annotations.NotNull;

public final class IclusionTestFactory {

    private IclusionTestFactory() {
    }

    @NotNull
    public static IclusionTrial trialWithTumors(@NotNull String acronym, @NotNull List<IclusionTumorLocation> tumorLocations,
            @NotNull List<IclusionTumorLocation> blacklist) {
        return trial(acronym, Lists.newArrayList(), tumorLocations, blacklist);
    }

    @NotNull
    public static IclusionTrial trial(@NotNull String acronym, @NotNull List<IclusionMutationCondition> mutationConditions) {
        return trial(acronym, mutationConditions, Lists.newArrayList(), Lists.newArrayList());
    }

    @NotNull
    public static IclusionTrial trialWithMutationsAndTumorLocation(@NotNull String acronym, @NotNull List<IclusionMutationCondition> mutationConditions,
            @NotNull List<IclusionTumorLocation> tumorLocations) {
        return trial(acronym, mutationConditions, tumorLocations, Lists.newArrayList());
    }

    @NotNull
    public static IclusionTrial trial(@NotNull String acronym, @NotNull List<IclusionMutationCondition> mutationConditions,
            @NotNull List<IclusionTumorLocation> tumorLocations, @NotNull List<IclusionTumorLocation> blacklist) {
        return ImmutableIclusionTrial.builder()
                .id("id")
                .acronym(acronym)
                .title("title")
                .eudra("eudra")
                .nct("nct")
                .ipn("ipn")
                .ccmo("ccmo")
                .mutationConditions(mutationConditions)
                .tumorLocations(tumorLocations)
                .blacklistedTumorLocations(blacklist)
                .build();
    }

    @NotNull
    public static IclusionMutationCondition or(@NotNull List<IclusionMutation> mutations) {
        return ImmutableIclusionMutationCondition.builder().logicType(IclusionMutationLogicType.OR).mutations(mutations).build();
    }

    @NotNull
    public static IclusionMutationCondition and(@NotNull List<IclusionMutation> mutations) {
        return ImmutableIclusionMutationCondition.builder().logicType(IclusionMutationLogicType.AND).mutations(mutations).build();
    }

    @NotNull
    public static IclusionMutation normalMutation() {
        return ImmutableIclusionMutation.builder().gene("gene").name("name").negation(false).build();
    }

    @NotNull
    public static IclusionMutation negatedMutation() {
        return ImmutableIclusionMutation.builder().gene("gene").name("name").negation(true).build();
    }
}