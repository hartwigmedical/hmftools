package com.hartwig.hmftools.serve.sources.actin;

import static org.junit.Assert.assertNotNull;

import java.io.IOException;

import com.google.common.io.Resources;

import org.junit.Test;

public class ActinReaderTest {

    private static final String EXAMPLE_TSV = Resources.getResource("actin/example.tsv").getPath();
    private static final String EXAMPLE_FILTER_TSV = Resources.getResource("actin/filter.tsv").getPath();

    @Test
    public void canReadFromExampleFiles() throws IOException {
        assertNotNull(ActinReader.read(EXAMPLE_TSV, EXAMPLE_FILTER_TSV));
    }
}