package com.hartwig.hmftools.orange.cohort.mapping;

import static org.junit.Assert.assertEquals;

import java.util.List;

import com.google.common.collect.ArrayListMultimap;
import com.google.common.collect.ListMultimap;
import com.google.common.collect.Lists;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.doid.DoidParents;
import com.hartwig.hmftools.common.doid.DoidTestFactory;
import com.hartwig.hmftools.orange.cohort.datamodel.ImmutableSampleData;
import com.hartwig.hmftools.orange.cohort.datamodel.SampleData;

import org.jetbrains.annotations.NotNull;
import org.junit.Test;

public class CohortMapperTest {

    @Test
    public void canMatchDoidsToCancerType() {
        CohortMapper mapper = createTestCohortMapper();

        assertEquals(CohortConstants.COHORT_OTHER, evaluate(mapper, "not a doid"));

        // DOID 1 maps to type 1 unless DOID = 1.2
        assertEquals("type 1", evaluate(mapper, "doid1.0"));
        assertEquals("type 1", evaluate(mapper, "doid1.1"));
        assertEquals(CohortConstants.COHORT_OTHER, evaluate(mapper, "doid1.2"));

        // DOID 2 maps to type 2 but only on exact match with DOID 2.1
        assertEquals(CohortConstants.COHORT_OTHER, evaluate(mapper, "doid2.0"));
        assertEquals("type 2", evaluate(mapper, "doid2.1"));
        assertEquals(CohortConstants.COHORT_OTHER, evaluate(mapper, "doid2.2"));

        // DOID 3/4/5 match to type 3/4/5 but DOID 5 has the lowest preference rank
        assertEquals("type 5", evaluate(mapper, "doid3", "doid4", "doid5"));
        // However, if only one DOID present, that should still win!
        assertEquals("type 4", evaluate(mapper, "doid4"));
        // DOID 4 and 7 have multiple mappings with same preference though!.
        assertEquals(CohortConstants.COHORT_OTHER, evaluate(mapper, "doid4", "doid7"));

        // DOID 6 maps to multiple cancer types on its own - wrong config.
        assertEquals(CohortConstants.COHORT_OTHER, evaluate(mapper, "doid6"));
    }

    @NotNull
    private static String evaluate(@NotNull CohortMapper mapper, @NotNull String... doids) {
        SampleData sample = ImmutableSampleData.builder().sampleId("TestSample").addDoids(doids).build();
        return mapper.cancerTypeForSample(sample);
    }

    private static CohortMapper createTestCohortMapper() {
        ListMultimap<String, String> relationship = ArrayListMultimap.create();
        relationship.put("doid1.2", "doid1.1");
        relationship.put("doid1.1", "doid1.0");

        relationship.put("doid2.2", "doid2.1");
        relationship.put("doid2.1", "doid2.0");
        DoidParents doidParentModel = DoidTestFactory.createDoidParents(relationship);

        ImmutableCohortMapping.Builder builder = ImmutableCohortMapping.builder().rule(MappingRule.DEFAULT).preferenceRank(3);

        List<CohortMapping> mappings = Lists.newArrayList();
        mappings.add(builder.cancerType("type 1").include(Sets.newHashSet("doid1.0")).exclude(Sets.newHashSet("doid1.2")).build());
        mappings.add(builder.cancerType("type 2").rule(MappingRule.EXACT_MATCH).include(Sets.newHashSet("doid2.1")).build());
        builder.rule(MappingRule.DEFAULT);

        mappings.add(builder.cancerType("type 3").include(Sets.newHashSet("doid3")).preferenceRank(2).build());
        mappings.add(builder.cancerType("type 4").include(Sets.newHashSet("doid4")).preferenceRank(3).build());
        mappings.add(builder.cancerType("type 6").include(Sets.newHashSet("doid7")).preferenceRank(3).build());
        mappings.add(builder.cancerType("type 5").include(Sets.newHashSet("doid5")).preferenceRank(1).build());

        mappings.add(builder.cancerType("type 6").include(Sets.newHashSet("doid6")).build());
        mappings.add(builder.cancerType("type 7").include(Sets.newHashSet("doid6")).build());

        return new CohortMapper(doidParentModel, mappings);
    }
}