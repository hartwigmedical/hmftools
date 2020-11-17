package com.hartwig.hmftools.serve.sources.iclusion.filter;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.iclusion.data.IclusionMutation;
import com.hartwig.hmftools.iclusion.data.IclusionMutationCondition;
import com.hartwig.hmftools.iclusion.data.IclusionMutationLogicType;
import com.hartwig.hmftools.iclusion.data.IclusionTrial;
import com.hartwig.hmftools.iclusion.data.ImmutableIclusionMutation;
import com.hartwig.hmftools.iclusion.data.ImmutableIclusionMutationCondition;
import com.hartwig.hmftools.iclusion.data.ImmutableIclusionTrial;

import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;
import org.junit.Test;

public class IclusionFilterTest {

    @Test
    public void canFilterIclusionTrials() {
        IclusionMutationCondition singleNormalOr = or(Lists.newArrayList(normalMutation()));
        IclusionMutationCondition singleNormalAnd = and(Lists.newArrayList(normalMutation()));
        IclusionMutationCondition singleNegatedOr = or(Lists.newArrayList(negatedMutation()));
        IclusionMutationCondition multipleNormalOr = or(Lists.newArrayList(normalMutation(), normalMutation()));
        IclusionMutationCondition multipleCombinedOr = or(Lists.newArrayList(normalMutation(), normalMutation(), negatedMutation()));

        IclusionTrial singleNormalOrTrial = trial("singleNormalOr", Lists.newArrayList(singleNormalOr));
        assertFalse(IclusionFilter.run(Lists.newArrayList(singleNormalOrTrial)).isEmpty());

        IclusionTrial singleNormalAndTrial = trial("singleNormalAnd", Lists.newArrayList(singleNormalAnd));
        assertTrue(IclusionFilter.run(Lists.newArrayList(singleNormalAndTrial)).isEmpty());

        IclusionTrial singleNegatedTrial = trial("singleNegatedOr", Lists.newArrayList(singleNegatedOr));
        assertTrue(IclusionFilter.run(Lists.newArrayList(singleNegatedTrial)).isEmpty());

        IclusionTrial multipleNormalOrTrial = trial("multipleNormalOr", Lists.newArrayList(multipleNormalOr));
        assertFalse(IclusionFilter.run(Lists.newArrayList(multipleNormalOrTrial)).isEmpty());

        IclusionTrial multipleCombinedOrTrial = trial("multipleCombinedOr", Lists.newArrayList(multipleCombinedOr));
        assertFalse(IclusionFilter.run(Lists.newArrayList(multipleCombinedOrTrial)).isEmpty());

        IclusionTrial noAcronymTrial = trial(Strings.EMPTY, Lists.newArrayList(singleNormalOr));
        assertTrue(IclusionFilter.run(Lists.newArrayList(noAcronymTrial)).isEmpty());

        // This trial should be kept since it has one multipleCombinedOr
        IclusionTrial complexTrial = trial("complex", Lists.newArrayList(singleNegatedOr, multipleCombinedOr));
        List<IclusionTrial> filteredComplex = IclusionFilter.run(Lists.newArrayList(complexTrial));
        assertEquals(1, filteredComplex.size());
        assertEquals(1, filteredComplex.get(0).mutationConditions().size());
    }

    @NotNull
    private static IclusionTrial trial(@NotNull String acronym, @NotNull List<IclusionMutationCondition> mutationConditions) {
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
    private static IclusionMutationCondition or(@NotNull List<IclusionMutation> mutations) {
        return ImmutableIclusionMutationCondition.builder().logicType(IclusionMutationLogicType.OR).mutations(mutations).build();
    }

    @NotNull
    private static IclusionMutationCondition and(@NotNull List<IclusionMutation> mutations) {
        return ImmutableIclusionMutationCondition.builder().logicType(IclusionMutationLogicType.AND).mutations(mutations).build();
    }

    @NotNull
    private static IclusionMutation normalMutation() {
        return ImmutableIclusionMutation.builder().gene("gene").name("name").negation(false).build();
    }

    @NotNull
    private static IclusionMutation negatedMutation() {
        return ImmutableIclusionMutation.builder().gene("gene").name("name").negation(true).build();
    }

}