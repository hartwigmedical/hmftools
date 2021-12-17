package com.hartwig.hmftools.serve.sources.actin.filter;

import static org.junit.Assert.assertEquals;

import java.io.IOException;
import java.util.List;

import com.google.common.io.Resources;

import org.junit.Test;

public class ActinFilterFileTest {

    private static final String TEST_ACTIN_FILTER_FILE = Resources.getResource("actin/filter.tsv").getPath();

    @Test
    public void canReadCkbFilterTsv() throws IOException {
        List<ActinFilterEntry> filterEntries = ActinFilterFile.read(TEST_ACTIN_FILTER_FILE);
        assertEquals(1, filterEntries.size());
    }

}