package com.hartwig.hmftools.serve.sources.actin.classification;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNotNull;
import static org.junit.Assert.assertTrue;

import java.util.Set;

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
    public void canExtractCorrectEvents() {
        ActinEntry activationAmplification =
                ActinTestFactory.builder().rule(ActinRule.ACTIVATION_OR_AMPLIFICATION_OF_GENE_X).gene("A").build();
        Set<String> events = ActinEventExtractor.extractEvents(activationAmplification);
        assertTrue(events.contains(ActinKeywords.ACTIVATION));
        assertTrue(events.contains(ActinKeywords.AMPLIFICATION));

        ActinEntry activation = ActinTestFactory.builder().rule(ActinRule.ACTIVATING_MUTATION_IN_GENE_X).gene("A").build();
        assertEquals(Sets.newHashSet(ActinKeywords.ACTIVATION), ActinEventExtractor.extractEvents(activation));

        ActinEntry promiscuousFusion = ActinTestFactory.builder().rule(ActinRule.FUSION_IN_GENE_X).gene("A").build();
        assertEquals(Sets.newHashSet(ActinKeywords.PROMISCUOUS_FUSION), ActinEventExtractor.extractEvents(promiscuousFusion));

        ActinEntry specificFusion = ActinTestFactory.builder().rule(ActinRule.SPECIFIC_FUSION_OF_X_TO_Y).mutation("EML4-ALK").build();
        assertEquals(Sets.newHashSet("EML4-ALK fusion"), ActinEventExtractor.extractEvents(specificFusion));

        ActinEntry inactivation = ActinTestFactory.builder().rule(ActinRule.INACTIVATION_OF_GENE_X).gene("A").build();
        assertEquals(Sets.newHashSet(ActinKeywords.INACTIVATION), ActinEventExtractor.extractEvents(inactivation));

        ActinEntry mutation = ActinTestFactory.builder().rule(ActinRule.MUTATION_IN_GENE_X_OF_TYPE_Y).gene("A").mutation("mut").build();
        assertEquals(Sets.newHashSet("mut"), ActinEventExtractor.extractEvents(mutation));

        ActinEntry amplification = ActinTestFactory.builder().rule(ActinRule.AMPLIFICATION_OF_GENE_X).gene("A").build();
        assertEquals(Sets.newHashSet(ActinKeywords.AMPLIFICATION), ActinEventExtractor.extractEvents(amplification));

        ActinEntry deletion = ActinTestFactory.builder().rule(ActinRule.DELETION_OF_GENE_X).gene("A").build();
        assertEquals(Sets.newHashSet(ActinKeywords.DELETION), ActinEventExtractor.extractEvents(deletion));

        ActinEntry wildtype = ActinTestFactory.builder().rule(ActinRule.WILDTYPE_OF_GENE_X).gene("A").build();
        assertEquals(Sets.newHashSet(ActinKeywords.WILDTYPE), ActinEventExtractor.extractEvents(wildtype));

        ActinEntry msi = ActinTestFactory.builder().rule(ActinRule.MSI_SIGNATURE).mutation("msi high").build();
        assertEquals(Sets.newHashSet("msi high"), ActinEventExtractor.extractEvents(msi));

        ActinEntry hrd = ActinTestFactory.builder().rule(ActinRule.HRD_SIGNATURE).mutation("HRD").build();
        assertEquals(Sets.newHashSet("HRD"), ActinEventExtractor.extractEvents(hrd));

        ActinEntry tmbHigh = ActinTestFactory.builder().rule(ActinRule.TMB_OF_AT_LEAST_X).mutation("TMB >= 10").build();
        assertEquals(Sets.newHashSet("TMB >= 10"), ActinEventExtractor.extractEvents(tmbHigh));

        ActinEntry tmlHigh = ActinTestFactory.builder().rule(ActinRule.TML_OF_AT_LEAST_X).mutation("TML >= 450").build();
        assertEquals(Sets.newHashSet("TML >= 450"), ActinEventExtractor.extractEvents(tmlHigh));

        ActinEntry tmlLow = ActinTestFactory.builder().rule(ActinRule.TML_OF_AT_MOST_X).mutation("TML < 290").build();
        assertEquals(Sets.newHashSet("TML < 290"), ActinEventExtractor.extractEvents(tmlLow));

        ActinEntry hla = ActinTestFactory.builder().rule(ActinRule.HAS_HLA_A_TYPE_X).mutation("HLA-02-1").build();
        assertEquals(Sets.newHashSet("HLA-02-1"), ActinEventExtractor.extractEvents(hla));
    }
}