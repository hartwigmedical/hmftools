package com.hartwig.hmftools.serve.sources.docm;

import static org.junit.Assert.assertNotNull;

import java.io.IOException;

import com.google.common.io.Resources;

import org.junit.Test;

public class DocmReaderTest {

    private static final String EXAMPLE_TSV = Resources.getResource("docm/example.tsv").getPath();

    @Test
    public void canReadDcomInputFile() throws IOException {
        assertNotNull(DocmReader.readAndCurate(EXAMPLE_TSV));
    }
}