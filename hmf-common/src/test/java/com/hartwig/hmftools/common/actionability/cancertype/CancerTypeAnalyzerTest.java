package com.hartwig.hmftools.common.actionability.cancertype;

import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import java.util.Map;
import java.util.Set;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.google.common.collect.Sets;

import org.junit.Test;

public class CancerTypeAnalyzerTest {

    @Test
    public void canMatchDOIDToPrimaryTumorLocation() {
        CancerTypeToDOIDMappingEntry cancerTypeToDOIDMappingEntry =
                ImmutableCancerTypeToDOIDMappingEntry.builder().cancerType("NSCLC").addDoids("0010").build();

        Map<String, Set<String>> primaryTumorLocationMappings = Maps.newHashMap();
        primaryTumorLocationMappings.put("Lung", Sets.newHashSet("0010"));

        CancerTypeAnalyzer cancerTypeAnalyzer = new CancerTypeAnalyzer(Lists.newArrayList(cancerTypeToDOIDMappingEntry),
                new PrimaryTumorToDOIDMapping(primaryTumorLocationMappings));

        assertFalse(cancerTypeAnalyzer.isCancerTypeMatch("Skin", "Lung"));
        assertFalse(cancerTypeAnalyzer.isCancerTypeMatch("Lung", "Lung"));
        assertTrue(cancerTypeAnalyzer.isCancerTypeMatch("NSCLC", "Lung"));
    }
}