package com.hartwig.hmftools.serve.sources.iclusion.filter;

import static com.hartwig.hmftools.serve.sources.iclusion.IclusionTestFactory.or;

import static org.junit.Assert.assertTrue;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.iclusion.datamodel.IclusionMutation;
import com.hartwig.hmftools.iclusion.datamodel.IclusionMutationCondition;
import com.hartwig.hmftools.iclusion.datamodel.IclusionTrial;
import com.hartwig.hmftools.iclusion.datamodel.IclusionTumorLocation;
import com.hartwig.hmftools.iclusion.datamodel.ImmutableIclusionMutation;
import com.hartwig.hmftools.iclusion.datamodel.ImmutableIclusionTrial;
import com.hartwig.hmftools.iclusion.datamodel.ImmutableIclusionTumorLocation;

import org.jetbrains.annotations.NotNull;
import org.junit.Test;

public class IclusionFilterTest {

    @Test
    public void canFilter() {
        List<IclusionFilterEntry> filterEntries = Lists.newArrayList();
        filterEntries.add(create(IclusionFilterType.FILTER_VARIANT_ON_GENE, "BRAF V600E"));
        filterEntries.add(create(IclusionFilterType.FILTER_EVENT_WITH_KEYWORD, "EXPRESSION"));

        IclusionFilter filter = new IclusionFilter(filterEntries);

        String location1 = "loc1";
        String loc1Doid1 = "loc1Doid1";
        String loc1Doid2 = "loc1Doid2";
        String location2 = "loc2";
        String loc2Doid1 = "loc2Doid2";

        IclusionTumorLocation loc1 =
                ImmutableIclusionTumorLocation.builder().primaryTumorLocation(location1).addDoids(loc1Doid1).addDoids(loc1Doid2).build();
        IclusionTumorLocation loc2 = ImmutableIclusionTumorLocation.builder().primaryTumorLocation(location2).addDoids(loc2Doid1).build();

        IclusionMutationCondition singleNormalOr = or(Lists.newArrayList((normalMutation())));

        IclusionTrial brafV600E =
                ImmutableIclusionTrial.builder()
                        .id("id")
                        .acronym("empty")
                        .title("title")
                        .eudra("eudra")
                        .nct("nct")
                        .ipn("ipn")
                        .ccmo("ccmo")
                        .mutationConditions(Lists.newArrayList(singleNormalOr))
                        .tumorLocations(Lists.newArrayList(loc1, loc2))
                        .build();

        assertTrue(filter.run(Lists.newArrayList(brafV600E)).isEmpty());

        IclusionTrial kras = ImmutableIclusionTrial.builder()
                .id("id")
                .acronym("empty")
                .title("title")
                .eudra("eudra")
                .nct("nct")
                .ipn("ipn")
                .ccmo("ccmo")
                .mutationConditions(Lists.newArrayList(singleNormalOr))
                .tumorLocations(Lists.newArrayList(loc1, loc2))
                .build();

        assertTrue(filter.run(Lists.newArrayList(kras)).isEmpty());

        filter.reportUnusedFilterEntries();
    }

    @NotNull
    public static IclusionMutation normalMutation() {
        return ImmutableIclusionMutation.builder().gene("gene").name("EXPRESSION").negation(false).build();
    }

    @NotNull
    private static IclusionFilterEntry create(@NotNull IclusionFilterType type, @NotNull String value) {
        return ImmutableIclusionFilterEntry.builder().type(type).value(value).build();
    }

}