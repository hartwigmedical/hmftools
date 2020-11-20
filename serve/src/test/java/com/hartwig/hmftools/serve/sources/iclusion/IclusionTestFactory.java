package com.hartwig.hmftools.serve.sources.iclusion;

import java.util.List;

import com.hartwig.hmftools.iclusion.data.IclusionMutation;
import com.hartwig.hmftools.iclusion.data.IclusionMutationCondition;
import com.hartwig.hmftools.iclusion.data.IclusionMutationLogicType;
import com.hartwig.hmftools.iclusion.data.IclusionTrial;
import com.hartwig.hmftools.iclusion.data.ImmutableIclusionMutation;
import com.hartwig.hmftools.iclusion.data.ImmutableIclusionMutationCondition;
import com.hartwig.hmftools.iclusion.data.ImmutableIclusionTrial;

import org.jetbrains.annotations.NotNull;

public final class IclusionTestFactory {

    private IclusionTestFactory() {
    }

    @NotNull
    public static IclusionTrial trial(@NotNull String acronym, @NotNull List<IclusionMutationCondition> mutationConditions) {
        return ImmutableIclusionTrial.builder()
                .id("id")
                .acronym(acronym)
                .title("title")
                .eudra("eudra")
                .nct("nct")
                .ipn("ipn")
                .ccmo("ccmo")
                .mutationConditions(mutationConditions)
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
