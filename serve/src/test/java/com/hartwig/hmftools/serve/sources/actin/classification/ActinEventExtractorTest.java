package com.hartwig.hmftools.serve.sources.actin.classification;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNotNull;

import com.google.common.collect.Sets;
import com.hartwig.hmftools.serve.sources.actin.ActinTestFactory;
import com.hartwig.hmftools.serve.sources.actin.reader.ActinEntry;
import com.hartwig.hmftools.serve.sources.actin.reader.ActinRule;

import org.junit.Test;

public class ActinEventExtractorTest {

    @Test
    public void canGenerateEventsForAllRules() {
        for (ActinRule rule : ActinRule.values()) {
            assertNotNull(ActinEventExtractor.extractEvents(ActinTestFactory.builder().rule(rule).build()));
        }
    }

    @Test
    public void canExtractEvent() {
        ActinEntry trial1 = ActinTestFactory.builder()
                .trial("trial1")
                .rule(ActinRule.MUTATION_IN_GENE_X_OF_TYPE_Y)
                .gene("A")
                .mutation("mut")
                .build();

        assertEquals(Sets.newHashSet("mut"), ActinEventExtractor.extractEvents(trial1));

        ActinEntry trial2 = ActinTestFactory.builder()
                .trial("trial1")
                .rule(ActinRule.WILDTYPE_OF_GENE_X)
                .gene("A")
                .mutation("wildtype")
                .build();

        assertEquals(Sets.newHashSet("wildtype"), ActinEventExtractor.extractEvents(trial2));

        ActinEntry trial3 = ActinTestFactory.builder()
                .trial("trial1")
                .rule(ActinRule.MSI_SIGNATURE)
                .gene("")
                .mutation("msi high")
                .build();

        assertEquals(Sets.newHashSet("MSI_high msi high"), ActinEventExtractor.extractEvents(trial3));

        ActinEntry trial4 = ActinTestFactory.builder()
                .trial("trial1")
                .rule(ActinRule.TML_OF_AT_LEAST_X)
                .gene("")
                .mutation("TML >= 450")
                .build();

        assertEquals(Sets.newHashSet("TML_high TML >= 450"), ActinEventExtractor.extractEvents(trial4));
    }
}