package com.hartwig.hmftools.common.serve.classification;

import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import com.google.common.collect.Lists;

import org.junit.Test;

public class CompositeEventMatcherTest {

    private static final String MATCHING_GENE = "match";
    private static final String NON_MATCHING_GENE = "no-match";

    @Test
    public void compositeEventClassifierWorks() {
        EventMatcher geneMatcher = (gene, event) -> gene.equals(MATCHING_GENE);
        EventMatcher noMatcher = (gene, event) -> false;

        EventMatcher noExclusions = new CompositeEventMatcher(Lists.newArrayList(), geneMatcher);
        assertTrue(noExclusions.matches(MATCHING_GENE, "any"));
        assertFalse(noExclusions.matches(NON_MATCHING_GENE, "any"));

        EventMatcher noExclusionsNoInclusions = new CompositeEventMatcher(Lists.newArrayList(), noMatcher);
        assertFalse(noExclusionsNoInclusions.matches(MATCHING_GENE, "any"));
        assertFalse(noExclusionsNoInclusions.matches(NON_MATCHING_GENE, "any"));

        EventMatcher excludingAndIncluding = new CompositeEventMatcher(Lists.newArrayList(geneMatcher), geneMatcher);
        assertFalse(excludingAndIncluding.matches(MATCHING_GENE, "any"));
        assertFalse(excludingAndIncluding.matches(NON_MATCHING_GENE, "any"));

        EventMatcher excludingAndNoInclusions = new CompositeEventMatcher(Lists.newArrayList(geneMatcher), noMatcher);
        assertFalse(excludingAndNoInclusions.matches(MATCHING_GENE, "any"));
        assertFalse(excludingAndNoInclusions.matches(NON_MATCHING_GENE, "any"));
    }
}