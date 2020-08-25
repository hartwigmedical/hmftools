package com.hartwig.hmftools.serve.dcom;

import static org.junit.Assert.assertEquals;

import java.io.IOException;
import java.util.List;

import com.google.common.io.Resources;

import org.junit.Test;

public class DcomFileReaderTest {

    private static final String EXAMPLE_TSV = Resources.getResource("dcom/example.tsv").getPath();

    @Test
    public void canReadDcomInputFile() throws IOException {
        List<DcomEntry> entries = DcomFileReader.readDcomFile(EXAMPLE_TSV);

        assertEquals(2, entries.size());

        assertEquals("MTOR", entries.get(0).gene());
        assertEquals("ENST00000361445", entries.get(0).transcript());
        assertEquals("p.R2505P", entries.get(0).proteinImpact());

        assertEquals("BTK", entries.get(1).gene());
        assertEquals("ENST00000372880", entries.get(1).transcript());
        assertEquals("p.E41K", entries.get(1).proteinImpact());
    }
}