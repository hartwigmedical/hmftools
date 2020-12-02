package com.hartwig.hmftools.serve.sources.iclusion.filter;

import static com.hartwig.hmftools.serve.sources.iclusion.IclusionTestFactory.and;
import static com.hartwig.hmftools.serve.sources.iclusion.IclusionTestFactory.negatedMutation;
import static com.hartwig.hmftools.serve.sources.iclusion.IclusionTestFactory.normalMutation;
import static com.hartwig.hmftools.serve.sources.iclusion.IclusionTestFactory.or;
import static com.hartwig.hmftools.serve.sources.iclusion.IclusionTestFactory.trial;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.iclusion.datamodel.IclusionMutationCondition;
import com.hartwig.hmftools.iclusion.datamodel.IclusionTrial;
import com.hartwig.hmftools.iclusion.datamodel.ImmutableIclusionMutation;

import org.apache.logging.log4j.util.Strings;
import org.junit.Test;

public class IclusionFilterTest {

    @Test
    public void canFilterIclusionTrials() {
        IclusionFilter filter = new IclusionFilter();
        IclusionMutationCondition singleNormalOr = or(Lists.newArrayList(normalMutation()));
        IclusionMutationCondition singleNormalAnd = and(Lists.newArrayList(normalMutation()));
        IclusionMutationCondition singleNegatedOr = or(Lists.newArrayList(negatedMutation()));
        IclusionMutationCondition multipleNormalOr = or(Lists.newArrayList(normalMutation(), normalMutation()));
        IclusionMutationCondition multipleCombinedOr = or(Lists.newArrayList(normalMutation(), normalMutation(), negatedMutation()));

        IclusionTrial singleNormalOrTrial = trial("singleNormalOr", Lists.newArrayList(singleNormalOr));
        assertFalse(filter.run(Lists.newArrayList(singleNormalOrTrial)).isEmpty());

        IclusionTrial singleNormalAndTrial = trial("singleNormalAnd", Lists.newArrayList(singleNormalAnd));
        assertTrue(filter.run(Lists.newArrayList(singleNormalAndTrial)).isEmpty());

        IclusionTrial singleNegatedTrial = trial("singleNegatedOr", Lists.newArrayList(singleNegatedOr));
        assertTrue(filter.run(Lists.newArrayList(singleNegatedTrial)).isEmpty());

        IclusionTrial multipleNormalOrTrial = trial("multipleNormalOr", Lists.newArrayList(multipleNormalOr));
        assertFalse(filter.run(Lists.newArrayList(multipleNormalOrTrial)).isEmpty());

        IclusionTrial multipleCombinedOrTrial = trial("multipleCombinedOr", Lists.newArrayList(multipleCombinedOr));
        assertFalse(filter.run(Lists.newArrayList(multipleCombinedOrTrial)).isEmpty());

        IclusionTrial noAcronymTrial = trial(Strings.EMPTY, Lists.newArrayList(singleNormalOr));
        assertTrue(filter.run(Lists.newArrayList(noAcronymTrial)).isEmpty());

        // This trial should be kept since it has one multipleCombinedOr
        IclusionTrial complexTrial = trial("complex", Lists.newArrayList(singleNegatedOr, multipleCombinedOr));
        List<IclusionTrial> filteredComplex = filter.run(Lists.newArrayList(complexTrial));
        assertEquals(1, filteredComplex.size());
        assertEquals(1, filteredComplex.get(0).mutationConditions().size());

        filter.reportUnusedFilterEntries();
    }

    @Test
    public void canFilterSpecificMutationsCorrectly() {
        IclusionFilter filter = new IclusionFilter();

        FilterKey firstKey = FilterFactory.MUTATION_KEYS_TO_FILTER.iterator().next();
        ImmutableIclusionMutation firstFilteredKey =
                ImmutableIclusionMutation.builder().gene(firstKey.gene()).name(firstKey.name()).negation(false).build();
        assertFalse(filter.include(firstFilteredKey));

        ImmutableIclusionMutation firstFilteredMutation = ImmutableIclusionMutation.builder()
                .gene("any")
                .name(FilterFactory.MUTATIONS_TO_FILTER.iterator().next())
                .negation(false)
                .build();

        assertFalse(filter.include(firstFilteredMutation));

        assertTrue(filter.include(ImmutableIclusionMutation.builder().gene("any").name("any").negation(false).build()));

        filter.reportUnusedFilterEntries();
    }
}