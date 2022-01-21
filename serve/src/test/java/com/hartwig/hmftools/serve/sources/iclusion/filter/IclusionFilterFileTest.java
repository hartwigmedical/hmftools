package com.hartwig.hmftools.serve.sources.iclusion.filter;

import static org.junit.Assert.assertEquals;

import java.io.IOException;
import java.util.List;

import com.google.common.io.Resources;

import org.junit.Test;

public class IclusionFilterFileTest {

    private static final String TEST_ICLUSION_FILTER_FILE = Resources.getResource("iclusion/filter.tsv").getPath();

    @Test
    public void canReadCkbFilterTsv() throws IOException {
        List<IclusionFilterEntry> filterEntries = IclusionFilterFile.read(TEST_ICLUSION_FILTER_FILE);
        assertEquals(1, filterEntries.size());
    }

}