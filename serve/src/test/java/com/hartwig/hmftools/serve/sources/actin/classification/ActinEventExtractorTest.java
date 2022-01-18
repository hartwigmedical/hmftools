package com.hartwig.hmftools.serve.sources.actin.classification;

import static org.junit.Assert.assertEquals;

import com.google.common.collect.Sets;
import com.hartwig.hmftools.serve.sources.actin.reader.ActinEntry;
import com.hartwig.hmftools.serve.sources.actin.reader.ActinRule;
import com.hartwig.hmftools.serve.sources.actin.reader.ImmutableActinEntry;

import org.junit.Test;

public class ActinEventExtractorTest {

    @Test
    public void canExtractEvent() {
        ActinEntry trial1 = ImmutableActinEntry.builder()
                .trial("trial1")
                .rule(ActinRule.MUTATION_IN_GENE_X_OF_TYPE_Y)
                .gene("A")
                .mutation("mut")
                .build();

        assertEquals(Sets.newHashSet("mut"), ActinEventExtractor.extractEvents(trial1));
    }
}