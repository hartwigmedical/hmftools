package com.hartwig.hmftools.serve.sources.actin.classification;

import static org.junit.Assert.assertEquals;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.serve.sources.actin.reader.ActinEntry;
import com.hartwig.hmftools.serve.sources.actin.reader.ActinRule;
import com.hartwig.hmftools.serve.sources.actin.reader.ImmutableActinEntry;

import org.junit.Test;

public class ActinEventAndGeneExtractorTest {

    @Test
    public void canExtractGene() {
        ActinEntry trial1 = ImmutableActinEntry.builder()
                .trial("trial1")
                .rule(ActinRule.ACTIVATING_FUSION_IN_GENE_X)
                .parameters(Lists.newArrayList("A")).build();
        assertEquals("A", ActinEventAndGeneExtractor.extractGene(trial1));
    }

    @Test
    public void canExtractEvent() {
        ActinEntry trial1 = ImmutableActinEntry.builder()
                .trial("trial1")
                .rule(ActinRule.MUTATION_IN_GENE_X_OF_TYPE_Y)
                .parameters(Lists.newArrayList("A", "mut")).build();
        assertEquals(Lists.newArrayList("mut"), ActinEventAndGeneExtractor.extractEvent(trial1));
    }
}