package com.hartwig.hmftools.serve.docm;

import static org.junit.Assert.assertEquals;

import java.io.IOException;
import java.util.List;

import com.google.common.io.Resources;

import org.junit.Test;

public class DocmFileReaderTest {

    private static final String EXAMPLE_TSV = Resources.getResource("docm/example.tsv").getPath();

    @Test
    public void canReadDcomInputFile() throws IOException {
        List<DocmEntry> entries = DocmFileReader.readDcomFile(EXAMPLE_TSV);

        assertEquals(2, entries.size());

        assertEquals("MTOR", entries.get(0).gene());
        assertEquals("ENST00000361445", entries.get(0).transcript());
        assertEquals("R2505P", entries.get(0).proteinAnnotation());

        assertEquals("BTK", entries.get(1).gene());
        assertEquals("ENST00000372880", entries.get(1).transcript());
        assertEquals("E41K", entries.get(1).proteinAnnotation());
    }
}