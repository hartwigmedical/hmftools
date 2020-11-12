package com.hartwig.hmftools.serve.sources.vicc.filter;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.serve.sources.vicc.ViccTestFactory;
import com.hartwig.hmftools.vicc.datamodel.Feature;
import com.hartwig.hmftools.vicc.datamodel.ImmutableFeature;
import com.hartwig.hmftools.vicc.datamodel.ViccEntry;
import com.hartwig.hmftools.vicc.datamodel.ViccSource;

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

        String nameToFilter = FilterFactory.FEATURE_KEYWORDS_TO_FILTER.iterator().next() + " filter me";
        Feature featureWithFilterKeyword = ImmutableFeature.builder().name(nameToFilter).build();
        Feature featureWithoutFilterKeyword = ImmutableFeature.builder().name("don't filter me").build();
        assertFalse(filter.include(ViccSource.CIVIC, featureWithFilterKeyword));
        assertTrue(filter.include(ViccSource.CIVIC, featureWithoutFilterKeyword));

        FilterKey keyToFilter = FilterFactory.FEATURE_KEYS_TO_FILTER.iterator().next();
        Feature featureToFilter = ImmutableFeature.builder().geneSymbol(keyToFilter.gene()).name(keyToFilter.name()).build();
        assertFalse(filter.include(keyToFilter.source(), featureToFilter));

        filter.reportUnusedFilterEntries();
    }
}