package com.hartwig.hmftools.serve.sources.actin.reader;

import static org.junit.Assert.assertEquals;

import java.io.IOException;
import java.util.List;

import com.google.common.collect.Lists;
import com.google.common.io.Resources;

import org.jetbrains.annotations.NotNull;
import org.junit.Test;

public class ActinFileReaderTest {

    private static final String EXAMPLE_TSV = Resources.getResource("actin/example.tsv").getPath();

    @Test
    public void canReadActinInputFile() throws IOException {
        List<ActinEntry> entries = ActinFileReader.read(EXAMPLE_TSV);

        assertEquals(3, entries.size());

        ActinEntry trial1 = findByTrial(entries, "Trial 1");
        assertEquals(ActinRule.AMPLIFICATION_OF_GENE_X, trial1.rule());
        assertEquals(Lists.newArrayList("CCND1"), trial1.parameters());


        ActinEntry trial2 = findByTrial(entries, "Trial 2");
        assertEquals(ActinRule.INACTIVATION_OF_GENE_X, trial2.rule());
        assertEquals(Lists.newArrayList("TP53"), trial2.parameters());

        ActinEntry trial3 = findByTrial(entries, "Trial 3");
        assertEquals(ActinRule.MUTATION_IN_GENE_X_OF_TYPE_Y, trial3.rule());
        assertEquals(Lists.newArrayList("CCND1", "3' UTR LOSS"), trial3.parameters());
    }

    @NotNull
    private static ActinEntry findByTrial(@NotNull List<ActinEntry> entries, @NotNull String trialToFind) {
        for (ActinEntry entry : entries) {
            if (entry.trial().equals(trialToFind)) {
                return entry;
            }
        }

        throw new IllegalStateException("Could not find actin entry for trial " + trialToFind);
    }
}