package com.hartwig.hmftools.serve.sources.actin.reader;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNull;

import java.io.IOException;
import java.util.List;

import com.google.common.io.Resources;

import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;
import org.junit.Test;

public class ActinFileReaderTest {

    private static final String EXAMPLE_TSV = Resources.getResource("actin/example.tsv").getPath();

    @Test
    public void canReadActinInputFile() throws IOException {
        List<ActinEntry> entries = ActinFileReader.read(EXAMPLE_TSV);

        assertEquals(3, entries.size());

        ActinEntry trial1 = findByTrial(entries, "Trial 1");
        assertEquals(ActinRule.TMB_OF_AT_LEAST_X, trial1.rule());
        assertEquals(Strings.EMPTY, trial1.gene());
        assertEquals("TMB >= 16", trial1.mutation());

        ActinEntry trial2 = findByTrial(entries, "Trial 2");
        assertEquals(ActinRule.INACTIVATION_OF_GENE_X, trial2.rule());
        assertEquals("TP53", trial2.gene());
        assertEquals(Strings.EMPTY, trial2.mutation());

        ActinEntry trial3 = findByTrial(entries, "Trial 3");
        assertEquals(ActinRule.MUTATION_IN_GENE_X_OF_TYPE_Y, trial3.rule());
        assertEquals("CCND1", trial3.gene());
        assertEquals("3' UTR LOSS", trial3.mutation());
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