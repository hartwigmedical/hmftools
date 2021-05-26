package com.hartwig.hmftools.serve.sources.ckb.filter;

import static org.junit.Assert.assertEquals;

import java.io.IOException;
import java.util.Set;

import com.google.common.io.Resources;

import org.junit.Test;

public class CkbFilterFileTest {

    private static final String TEST_CKB_FILTER_FILE = Resources.getResource("ckb_filter/ckb_filters.tsv").getPath();

    @Test
    public void canReadCkbFilterTsv() throws IOException {
        Set<CkbFilterEntry> filterEntries = CkbFilterFile.read(TEST_CKB_FILTER_FILE);
        assertEquals(2, filterEntries.size());
    }
}