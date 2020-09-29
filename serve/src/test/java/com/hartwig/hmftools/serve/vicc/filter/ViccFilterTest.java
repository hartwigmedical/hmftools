package com.hartwig.hmftools.serve.vicc.filter;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.serve.vicc.ViccTestFactory;
import com.hartwig.hmftools.vicc.datamodel.Feature;
import com.hartwig.hmftools.vicc.datamodel.ImmutableFeature;
import com.hartwig.hmftools.vicc.datamodel.ViccEntry;

import org.junit.Test;

public class ViccFilterTest {

    @Test
    public void canFilterOncogenicEvents() {
        ViccEntry oncogenic = ViccTestFactory.testViccEntryWithOncogenic("Oncogenic");
        ViccEntry benign = ViccTestFactory.testViccEntryWithOncogenic("Inconclusive");

        ViccFilter filter = new ViccFilter();
        List<ViccEntry> filteredEntries = filter.run(Lists.newArrayList(oncogenic, benign));
        assertEquals(1, filteredEntries.size());
        assertTrue(filteredEntries.contains(oncogenic));

        filter.reportUnusedFilterEntries();
    }

    @Test
    public void canFilterIndividualFeatures() {
        ViccFilter filter = new ViccFilter();

        String featureToFilter = FilterFactory.FEATURE_KEYWORDS_TO_FILTER.iterator().next() + " filter me";
        Feature featureWithFilterKeyword = ImmutableFeature.builder().name(featureToFilter).build();
        Feature featureWithoutFilterKeyword = ImmutableFeature.builder().name("don't filter me").build();
        assertFalse(filter.include(featureWithFilterKeyword));
        assertTrue(filter.include(featureWithoutFilterKeyword));
    }
}