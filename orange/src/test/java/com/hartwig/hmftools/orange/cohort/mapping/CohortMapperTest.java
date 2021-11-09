package com.hartwig.hmftools.orange.cohort.mapping;

import static org.junit.Assert.assertNotNull;

import java.util.List;

import com.google.common.collect.ArrayListMultimap;
import com.google.common.collect.ListMultimap;
import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.doid.DoidParents;
import com.hartwig.hmftools.common.doid.DoidTestFactory;

import org.junit.Test;

public class CohortMapperTest {

    @Test
    public void canMatchDoidsToCancerType() {
        ListMultimap<String, String> relationship = ArrayListMultimap.create();
        relationship.put("chain1", "chain1.1");
        relationship.put("chain1.1", "chain1.2");

        relationship.put("chain2", "chain2.1");
        relationship.put("chain2.1", "chain2.2");
        DoidParents doidParentModel = DoidTestFactory.createDoidParents(relationship);

        List<CohortMapping> mappings = Lists.newArrayList();
        mappings.add(ImmutableCohortMapping.builder()
                .cancerType("type 1")
                .mappingRule(MappingRule.DEFAULT)
                .addInclude("chain1")
                .addExclude("chain2")
                .preferenceRank(1)
                .build());

        CohortMapper mapper = new CohortMapper(doidParentModel, mappings);
        assertNotNull(mapper);
    }
}