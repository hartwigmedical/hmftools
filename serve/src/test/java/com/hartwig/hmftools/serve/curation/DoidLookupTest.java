package com.hartwig.hmftools.serve.curation;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNull;

import java.util.Map;
import java.util.Set;

import com.google.common.collect.Maps;
import com.google.common.collect.Sets;

import org.junit.Test;

public class DoidLookupTest {

    @Test
    public void canLookupDoids() {
        String testCancerType = "cancerA";
        Set<String> testDoids = Sets.newHashSet("123");

        Map<String, Set<String>> testMapping = Maps.newHashMap();
        testMapping.put(testCancerType, testDoids);

        DoidLookup doidLookup = new DoidLookup(testMapping);
        assertEquals(testDoids, doidLookup.lookupDoidsForCancerType(testCancerType));
        assertNull(doidLookup.lookupDoidsForCancerType("does not exist"));

        doidLookup.evaluate();
    }
}