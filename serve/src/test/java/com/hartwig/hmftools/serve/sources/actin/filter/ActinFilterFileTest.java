package com.hartwig.hmftools.serve.sources.actin.filter;

import static org.junit.Assert.assertEquals;

import java.io.IOException;
import java.util.List;

import com.google.common.io.Resources;

import org.junit.Test;

public class ActinFilterFileTest {

    private static final String EXAMPLE_FILTER_TSV = Resources.getResource("actin/filter.tsv").getPath();

    @Test
    public void canReadActinFilterTsv() throws IOException {
        List<ActinFilterEntry> filterEntries = ActinFilterFile.read(EXAMPLE_FILTER_TSV);

        assertEquals(1, filterEntries.size());
    }
}